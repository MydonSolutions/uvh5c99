using YAML
using Blio
using Dates
using XGPU
using HDF5
using ERFA
using Geodesy
using EarthOrientation
using RadioInterferometry

using Rawx


"""
!
! Copied from Rawx, in order to avoid its sorts which aren't done in UVH5C99:
!  - `pols = sort(unique(getindex.(input_map, 2)))`
!  - `antidxs = sort(unique(getindex.(input_map, 1)))`
!


    (pol_array, bl_array, inpair_maps) = generateInputPairMap(antennas, input_map)

For the given `antennas` and `input_map`, create polarization array
`pol_array::Array{Int,1}` containing UVH5/AIPS cross polarization values, a
baseline array `bl_array::Array{Tuple{Int,Int},1}` containing all pairs of all
antenna found in `input_map`, and an input pair mapping list
`inpair_maps::Array{InputPairMapping,1}`.
"""
function generateInputPairMap(antennas::AbstractArray{<:AbstractDict{Symbol,Any},1},
                              input_map::AbstractArray{<:AbstractArray{String,1},1})
  # Make sure we do not have more inputs than supported by xGPU
  num_inputs = length(input_map)
  @assert num_inputs <= xgpuinfo.nstation * xgpuinfo.npol "input_map too long for xGPU"

  # Remap antennas in input_map from name (String) or number (Int) to index
  # into antennas array and map polarizations to lowercase Symbols
  input_map = map(input_map) do (ant, pol)
    Rawx.mapantpol(ant, pol, antennas)
  end

  # Determine unique polarizations in input_map.
  #! pols = sort(unique(getindex.(input_map, 2)))
  pols = unique(getindex.(input_map, 2))
  npols = length(pols)
  @assert (npols == 1 || npols == 2) "npols must be 1 or 2 (not $npols)"

  # pol_array contains UVH5/AIPS polarization values for all possible
  # polarization pairings.  This is the order in which the polarizations will
  # be stored for a given antenna pair.
  pol_array = [
    Rawx.POLMAP[p1, p2] for p1 in pols for p2 in pols
  ]

  # pol_map (not to be confused with POLMAP!) is a Dict mapping (:p1,:p2)
  # polarization tuples to "polidx" polarization indexes (i.e. indexes into
  # polarization_array).
  pol_map = Dict(map(i->(Rawx.POLMAP[pol_array[i]]=>i), eachindex(pol_array)))

  # GUPPI RAW considers each antenna to be dual polarization.  xGPU can
  # be compiled for single or dual polarization inputs.  UVH5 files can
  # be single-polarization or dual-polarization.  In reality, the
  # actual antenna-polarizations going to each input can be a completely
  # arbitrary combination and arrangement of antennas and polarizations.
  #
  # We need to build a list of all possible antenna pairings (aka baselines)
  # given the set of unique antennas that appear in input_map.  It is possible
  # that some baselines will not have all cross polarizations present because
  # some antennas may not have both polarizations present.  The ordering of
  # antennas comprising a baseline is also important.  Because the
  # antenna-polarizations may be connected to arbitrary inputs of the
  # correlator, it is possible that some cross-polarizations will be conjugated
  # differently from other cross-polarizations of the same baseline.  Since
  # there is only one UVW vector for the baseline, these differences must be
  # resolved by conjugating some of the cross-polarizations of these
  # "conflicted" baselines to ensure that all cross-polarizations are
  # conjugated consistently.
  #
  # xGPU computes conj(in1)*in2 where in1 <= in2.  By convention, we constrain
  # baselines (antidx1, antidx2) such that antidx1 <= antidx2.  Note that
  # antidx is the index into the `antennas` array for a given antenna.  This
  # approach relies on the ordering of the antennas array rather than the
  # actual antenna names or numbers, which are essentially opaque labels with
  # no inherent meaning.
  #
  # UVH5 follows the conjugation convention of the Radio Interferometry
  # Measurement Equation (RIME) which defines the visibility product of
  # baseline (antidx1, antidx2) to be ``z_{antidx1}*conj(z_{antidx2})``.
  #
  # If xGPU input pair (in1, in2) with in1 < in2 corresponds to antenna index
  # pair (antidx1, antidx2) with antidx1 <= antidx2, then we consider the
  # visibility product from xGPU to be in need of conjugation before being
  # written to the UVH5 dataset. If we find an xGPU input pair that corresponds
  # to antenna index pair (antidx1, antidx2) with antidx1 > antidx2, then we
  # consider the visibility product be properly conjugated.
  #
  # Auto correlations for dual-polarization systems are handled as a special
  # case, because we write out the extra redundant cross-polarization product.
  # Extra care is required to ensure that these cross-polarization products are
  # properly conjugated.
  #
  # The ordering of the baselines in the dataset is autocorrelations followed
  # by cross correlation.  First we build a baseline array, bl_array, that will
  # be the source of the ant_1_array and ant_2_array contents.  Then we use
  # this list to build a baseline map, `bl_map` that maps (antidx1, antidx2)
  # baselines to a sequential baseline index.  This is used to look up the
  # value for a `bl_idx` field for a given antidx pair (for a given xGPU input
  # pair) when generating the list of InputPairMappings.

  #! antidxs = sort(unique(getindex.(input_map, 1)))
  antidxs = unique(getindex.(input_map, 1))
  nants = length(antidxs)
  nbls = ((nants + 1) * nants) ÷ 2

  bl_array = Array{Tuple{Int,Int},1}(undef, nbls)
  i = 0
  # autos
  for ai in antidxs
    i += 1
    bl_array[i] = (ai, ai)
  end
  # crosses
  for aii1 in 1:nants-1
    for aii2 in aii1+1:nants
      i += 1
      bl_array[i] = (antidxs[aii1], antidxs[aii2])
    end
  end
  # bl_map maps (antidx1, antodx2) tuples to baseline index
  bl_map = Dict(map(i->(bl_array[i]=>i), eachindex(bl_array)))

  # Now we can make list of input pair mappings.
  num_inpairs = (num_inputs+1) * num_inputs ÷ 2
  inpair_maps = Array{Rawx.InputPairMapping,1}(undef, num_inpairs)
  i = 0
  for in1 in 1:num_inputs
    ant1, pol1 = input_map[in1]
    for in2 in in1:num_inputs
      ant2, pol2 = input_map[in2]
      isauto = (ant1 == ant2)

      # TODO Need to check pol order as well (at least for autos)
      if haskey(bl_map, (ant1, ant2))
        blidx = bl_map[ant1, ant2]
        polidx = pol_map[pol1, pol2]
        needsconj = true
      else
        blidx = bl_map[ant2, ant1]
        polidx = pol_map[pol2, pol1]
        needsconj = false
      end

      i += 1
      inpair_maps[i] = Rawx.InputPairMapping((XGPU.xgpuInputPairIndex(in1, in2),
                                         blidx, polidx, isauto, needsconj))
    end # for in2
  end # for in1

  (pol_array, bl_array, inpair_maps)
end


"""
    uvh5debug_fromraw(rawstem::AbstractString,
              uvh5name::AbstractString;
              telinfo::Union{AbstractString,AbstractDict{Symbol,Any}},
              obsinfo::Union{AbstractString,AbstractDict{Symbol,Any}},
              subints_per_dump::Integer=0,
              inttime_per_dump::Real=0
             )::Nothing

Feeds data from GUPPI RAW files corresponding to `rawstem` to xGPU and outputs
visibilities in UVH5 format to file `uvh5name`.  If `rawstem` has a `.NNNN.raw`
suffix, it will be removed.

Metadata for the telescope array is provided in `telinfo`.  Metadata for the
observation that produced the RAW files is provided in `obsinfo`.  These
can be passed as `AbstractDict{Symbol,Any}`, in which case they are used as is,
or as `AbstractString` which are treated as the name of a YAML file from which
the metadata will be loaded.

Integration time is specified by `subints_per_dump` which is the number of xGPU
input buffers to correlate and integrate before dumping, a value of 0 means
calculate from `inttime_per_dump`, values less than zero mean dump every
`-subints_per_dump` blocks.  If both `subints_per_dump` and `inttime_per_dump`
are zero, then dump once per block.  `inttime_per_dump` must be given in
seconds.

A "subint" or "sub-integration" is the integration time granularity set by
XGPU's `xgpuinfo.ntime` value.  The number of time samples per GUPPI RAW block
must be divisible by `xgpuinfo.ntime`.  The integration time of a
sub-integration is given by `TBIN * xgpuinfo.ntime`.

The `telinfo` and `obsinfo` dictionaries are ultimately merged together, so
their separation here is really just for the caller's convenience.  Additional
metadata from the GUPPI RAW header is also added, but not if the caller has
provided explicit values to override the GUPPI RAW header (e.g. corrections to
the source name, sky coordinates, etc.).

Still TODO: Document `telinfo` and `obsinfo` more.
"""
function uvh5debug_fromraw(rawstem::AbstractString,
                   uvh5name::AbstractString;
                   telinfo::Union{AbstractString,AbstractDict{Symbol,Any}},
                   obsinfo::Union{AbstractString,AbstractDict{Symbol,Any}},
                   subints_per_dump::Integer=0,
                   inttime_per_dump::Real=0
                  )::Nothing

  @assert !ispath(uvh5name) "output file $uvh5name exists"
  @assert inttime_per_dump >= 0 "inttime_per_dump cannot be negative"

  # Keyword args don't participate in multiple dispatch so we have to handle the
  # AbstractString case here ourselves.
  if telinfo isa AbstractString
    telinfo = YAML.load_file(telinfo, dicttype=Dict{Symbol,Any})
  end
  if obsinfo isa AbstractString
    obsinfo = YAML.load_file(obsinfo, dicttype=Dict{Symbol,Any})
  end

  # Merge telinfo and obsinfo into metadata
  metadata = merge(telinfo, obsinfo)

  pol_array, bl_array, inpair_maps = generateInputPairMap(metadata[:antennas],
                                                          metadata[:input_map])
  # @info "inpair_maps length $idxs_length"
  # println(inpair_maps)
  # exit(1)

  # GUPPI RAW considers each antenna to be dual polarization.  xGPU can
  # be compiled for single or dual polarization inputs.  UVH5 files can
  # be single polarization of dual-polarization.  In reality, the
  # actual antenna-polarizations going to each input can be a completely
  # arbitrary combination and arrangement of antennas and
  # polarizations.  We must determine nants_data from the number of
  # unique antennas that appear in input_map.
  unique_ants = sort(unique(getindex.(metadata[:input_map], 1)))
  nants_data = length(unique_ants)
  nbls = ((nants_data + 1) * nants_data) ÷ 2

  # Set metadata polarization_array
  metadata[:polarization_array] = pol_array
  npols = length(pol_array)

  # Create Arrays for extendable datasets that do NOT depend on RAW headers or
  # data.

  # The a1s and a2s arrays are antenna_numbers, NOT indices (despite the
  # wording in the 2018-11-28 version of the UVH5 memo).  The antenna indexes
  # in bl_array are 1-based indexes into antenna_numbers and antenna_names.
  antenna_numbers = Int.(getindex.(metadata[:antennas], :number))
  a1s = antenna_numbers[getindex.(bl_array, 1)]
  a2s = antenna_numbers[getindex.(bl_array, 2)]
  antpos_enu = Rawx.getantpos(metadata, :enu)
  # Setup variables for ENU->UVW conversion
  lat_rad = Rawx.dms2rad(metadata[:latitude])
  lon_rad = Rawx.dms2rad(metadata[:longitude])
  alt_m = metadata[:altitude]
  antpos_uvw = similar(antpos_enu)
  bl_uvws = Array{Float64,2}(undef, 3, nbls)

  # Trim off .NNNN.raw suffix as a convenience
  rawstem = replace(rawstem, r"\.\d+\.raw$"=>"")
  # Take hash of rawstem for use as log message ID
  rawstemid = hash(rawstem)

  @info "working stem: $rawstem"

  rawhdr::GuppiRaw.Header = GuppiRaw.Header()
  uvh5 = nothing
  uvh5hdr = nothing
  uvh5data = nothing
  uvh5debug = nothing
  xgpuctx = nothing
  array_h = nothing
  matrix_h = nothing
  rawdata = nothing
  tau = 0.0
  timejd = 0.0
  hob = 0.0
  dob = 0.0
  lst = 0.0
  visdata = nothing
  nblts = 0
  ntimes = 0
  nfreqs = 0

  ra_rad = 0.0
  dec_rad = 0.0
  recalc_uvws = false

  rawseq = 0
  subints_sent = 0 # Total number of subints sent to XGPU
  subints_zero = 0 # Number of subints of zeros sent to XGPU
  subints_per_block = 0

  # PKTIDX of next block to process
  pktidx = -1
  # PKTIDX per block
  piperblk = 0

  rawname = Rawx.mkrawname(rawstem, rawseq)
  if !isfile(rawname)
    @error "opening file: $rawname [not found]"
    return
  end

  @info "opening file: $rawname"
  open(rawname) do rawio

    # While header
    while read!(rawio, rawhdr)
      # Get PKTIDX from header
      hdr_pktidx = rawhdr[:pktidx]::Integer

      # If no pktidx, this is the first block of first file
      if pktidx == -1
        # Get PKTSTART (use PKTIDX if not found)
        hdr_pktstart = get(rawhdr, :pktstart, hdr_pktidx)
        piperblk = get(rawhdr, :piperblk, 0)
        piperblk == 0 || @info "found piperblk = $piperblk"
        # Check timeliness of first block
        if hdr_pktidx == hdr_pktstart
          # All good!
          pktidx = hdr_pktidx
        elseif hdr_pktidx > hdr_pktstart
          if piperblk != 0 && (hdr_pktidx - hdr_pktstart) % piperblk == 0
            # Missed integer blocks at start, use PKTSTART as pktidx
            @warn "missed integer blocks at start" _module=nothing _file=nothing
            pktidx = hdr_pktstart
          else
            if piperblk != 0
              @warn "missed non-integer blocks at start" _module=nothing _file=nothing
            else
              @warn "data starts after pktstart" _module=nothing _file=nothing
            end
            pktidx = hdr_pktidx
          end
        else # hdr_pktidx < hdr_pktstart
          @warn "data starts before pktstart" _module=nothing _file=nothing
          pktidx = hdr_pktidx
        end
      # elseif no piperblk, calculate it from this block's pktidx and previous pktidx
      elseif piperblk == 0
        piperblk = hdr_pktidx - pktidx
        pktidx = hdr_pktidx
        @info "calculated piperblk = $piperblk"
      else
        # We have pktidx and piperblk.
        # Here are the possible situations:
        if hdr_pktidx > pktidx
          # Missing blocks.  Data for missing blocks will be handled below,
          # but show error message here if gap is not a multiple of piperblk
          # and increment pktidx to the smallest value is a multiple of
          # piperblk from hdr_pktidx.
          if (hdr_pktidx - pktidx) % piperblk != 0
            @error("non-integer number of missing blocks detected", rawstem,
                    maxlog=1, _id=Symbol(:nonint_missing_, rawstemid),
                    _module=nothing, _file=nothing)
            # "Catch up" to hdr_pktidx "phase"
            pktidx += (hdr_pktidx - pktidx) % piperblk
          end
        elseif pktidx > hdr_pktidx
          # Show error message
          @error("pktidx incremented by less than piperblk",
                  maxlog=1, _id=Symbol(:small_dpktidx_, rawstemid),
                  _module=nothing, _file=nothing)
          # If header does not have PIPERBLK, then earlier piperblk
          # calculation was wrong (as can happen if blocks were missed between
          # very first two blocks)
          if !haskey(rawhdr, :piperblk)
            # Adjust PIPERBLK by "overshoot"
            piperblk -= pktidx - hdr_pktidx
            @info "calculated new piperblk = $piperblk"
          end
          # Use hdr_pktidx as pktidx
          pktidx = hdr_pktidx
        end
      end

      # Setup UVH5 file if not already setup
      if isnothing(uvh5)
        # Validate that RAW file is compatible with xGPU params
        ntime = Blio.ntime(rawhdr)
        # Make sure that RAW ntime is divisible by xGPU NTIME
        @assert ntime % xgpuinfo.ntime == 0 "RAW ntime not divisible by xgpuinfo.ntime"
        # We currently require that RAW nants <= xGPU NSTATION
        @assert rawhdr.nants <= xgpuinfo.nstation
        # We currently require that all xGPU inputs are mappable to antpols
        @assert xgpuinfo.nstation * xgpuinfo.npol >= length(metadata[:input_map])

        subints_per_block = ntime ÷ xgpuinfo.ntime
        inttime_per_subint = rawhdr.tbin * xgpuinfo.ntime
        if subints_per_dump == 0
          if inttime_per_subint == 0
            # Default mode of one dump per GUPPI RAW block
            subints_per_dump = subints_per_block
          else
            # Compute subints_per_dump from inttime_per_subint
            subints_per_dump = Int(round(inttime_per_dump/inttime_per_subint))
          end
        elseif subints_per_dump < 0
          subints_per_dump = subints_per_block * -subints_per_dump
        end

        # Calculate actual inttime_per_dump
        inttime_per_dump = inttime_per_subint * subints_per_dump
        @info "integration details" subints_per_block subints_per_dump (
                                    inttime_per_subint) inttime_per_dump

        # Add to metadata based on GUPPI RAW headers unless the caller has
        # provided overrides.

        # Use GUPPI RAW `src_name` header as `object_name` unless the user
        # provided `object_name` in metadata.
        get!(metadata, :object_name, get(rawhdr, :src_name, "UNKNOWN"))
        # TODO Provide more convenient way to allow for freq_array override.
        get!(metadata, :freq_array, collect(chanfreqs(rawhdr)*1e6))
        nfreqs = length(metadata[:freq_array])
        # We currently require that nfreqs <= xGPU nfrequency
        @assert nfreqs <= xgpuinfo.nfrequency "RAW nfreqs ($nfreqs) > XGPU nfrequency ($(xgpuinfo.nfrequency))"

        get!(metadata, :channel_width, abs(get(rawhdr, :chan_bw, 0)) * 1e6)
        if !isa(metadata[:channel_width], AbstractVector)
          metadata[:channel_width] = fill(metadata[:channel_width], nfreqs)
        end

        # If RAW file has RA and DEC values in the header, phase_type is
        # "phased", otherwise it is "drift".
        if haskey(rawhdr, :ra) && haskey(rawhdr, :dec)
          metadata[:phase_type] = "phased"
          ra_rad = get!(metadata, :phase_center_ra, deg2rad(rawhdr.ra))
          dec_rad = get!(metadata, :phase_center_dec, deg2rad(rawhdr.dec))
          # TODO Handle case when phase_center_epoch is not 2000.0
          get!(metadata, :phase_center_epoch, 2000.0)
          get!(metadata, :phase_center_frame, "icrs")
          recalc_uvws = true
        else
          # Drift scan is phased to zenith, so antpos_uvw is identical to antpos_enu
          metadata[:phase_type] = "drift"
          Rawx.calc_uvws!(bl_uvws, antpos_enu, bl_array)
          recalc_uvws = false
        end

        # Create and initialize uvh5 file
        uvh5 = Rawx.uvh5create(uvh5name; metadata...)
        uvh5hdr = uvh5["Header"]
        uvh5data = uvh5["Data"]

        uvh5debug = create_group(uvh5, "Debug")
        create_dataset(uvh5debug, "dut1", Float64, ((0,),(-1,)),
                chunk=(1,))
        create_dataset(uvh5debug, "ra_rad", Float64, ((0,),(-1,)),
                chunk=(1,))
        create_dataset(uvh5debug, "dec_rad", Float64, ((0,),(-1,)),
                chunk=(1,))

        # Init xGPU and get context
        xgpuctx=xgpuInit()
        array_h=xgpuInputArray(xgpuctx);
        matrix_h=xgpuOutputArray(xgpuctx);

        matrix_h_size = size(matrix_h)
        @info "matrix_h size $matrix_h_size" 
        create_dataset(uvh5debug, "matrix_h", Complex{Int32},
                ((matrix_h_size..., 0), (matrix_h_size..., -1)),
                chunk=(matrix_h_size..., 1),
              )

        # Create rawdata Array
        rawdata = Array(rawhdr)

        # Initialize time variables
        tau = Float64(inttime_per_dump)
        # Initialize timejd to half of tau *before* start time of file because
        # each dump will increment times by tau.
        timejd = Rawx.starttime(rawhdr, pktidx) - tau/2 / ERFA.DAYSEC
        # Create visdata Array
        visdata = zeros(eltype(matrix_h), npols, nfreqs, nbls)
      end # if isnothing(uvh5) # lazy init

      # Avoid infinite loop
      if pktidx != hdr_pktidx
        if piperblk <= 0 || (hdr_pktidx - pktidx) % piperblk != 0
          @error("pktidx snafu", rawstem, pktidx, hdr_pktidx, piperblk,
                  maxlog=1, _id=Symbol(:pktidx_snafu_, rawstemid),
                  _module=nothing, _file=nothing)
          pktidx = hdr_pktidx
        end
      end

      # Loop from pktidx to hdr_pktidx, which should be a single iteration,
      # but could be more if blocks were missed.
      while pktidx <= hdr_pktidx # block loop
        if pktidx < hdr_pktidx
          @warn("missing blocks detected", rawstem, maxlog=1,
                _id=Symbol(:missing_blocks_, rawstemid),
                _module=nothing, _file=nothing)
          # Send block of zeros
          fill(rawdata, 0)
          subints_zero += subints_per_block
        else
          # Try to read block
          # On error, print message, set data_valid=false and break out of loop
          try
            read!(rawio, rawdata)
          catch ex
            @warn "got $ex while reading $rawname [done]" _module=nothing _file=nothing
            data_valid = false
            break # out of block loop
          end
        end
  #=
        # Debug, zero out data, then set test values
        rawdata .= 0
        rawdata[1,1,1,1] = 1+2im # indexed as (pol, time, chan, ant)
        rawdata[2,1,1,1] = 3+4im # indexed as (pol, time, chan, ant)
        rawdata[2,1,1,2] = 3+4im # indexed as (pol, time, chan, ant)
  =#

        # For subint in block
        for i in 1:subints_per_block
          # Swizzle subint data into XGPU input buffer array_h
          rawSwizzleInput!(array_h, rawdata, xgpuinfo.ntime*(i-1)+1, xgpuinfo,
                          rawhdr.nants, nfreqs)
          # send subint to GPU
          xgpuCudaXengine(xgpuctx, XGPU.SYNCOP_SYNC_TRANSFER)
          # Count subints
          subints_sent += 1

          # If nsent == subints_per_dump
          if subints_sent % subints_per_dump == 0
            # Dump integration buffer to matrix_h
            xgpuDumpDeviceIntegrationBuffer(xgpuctx)
            # Clear integration buffer on device
            xgpuClearDeviceIntegrationBuffer(xgpuctx)

            # Increment timejd
            timejd += tau / ERFA.DAYSEC

            # Get interpolated ut1utc value
            ut1utcsec = getΔUT1(timejd)
            ut1utcday = ut1utcsec / ERFA.DAYSEC

            # Calculate LST.  This uses UT1 for both UT1 and TT, which results
            # in an error on the order of 100 microarcseconds or approximately
            # 7 microseconds.
            lst = ERFA.gst06a(timejd, ut1utcday, timejd, ut1utcday) + lon_rad

            # Update uvws if we have the required info
            if recalc_uvws
              aob, zob, hob, dob, rob, eo = ERFA.atco13(
                ra_rad, dec_rad,  # right ascension, declination
                0, 0, 0, 0,       # proper motion ra/dec, parallax, rv,
                timejd, 0,        # jd1 (UTC), jd2 (UTC)
                ut1utcsec,        # DUT1 (UT1-UTC) at timejd
                lon_rad, lat_rad, alt_m, # observer location
                0, 0,             # polar motion (UVH5 does not record this)
                0, 0, 0, 0)       # atmosphereic refraction

              enu2uvw!(antpos_uvw, antpos_enu, hob, dob, lat_rad)

              Rawx.calc_uvws!(bl_uvws, antpos_uvw, bl_array)
            end

            # Extend datasets
            Rawx.extend(uvh5hdr["ant_1_array"])
            Rawx.extend(uvh5hdr["ant_2_array"])
            Rawx.extend(uvh5hdr["uvw_array"])
            Rawx.extend(uvh5hdr["time_array"])
            Rawx.extend(uvh5hdr["lst_array"])
            Rawx.extend(uvh5hdr["integration_time"])
            Rawx.extend(uvh5data["flags"])
            Rawx.extend(uvh5data["nsamples"])
            Rawx.extend(uvh5data["visdata"])
            
            Rawx.extend(uvh5debug["dut1"])
            Rawx.extend(uvh5debug["ra_rad"])
            Rawx.extend(uvh5debug["dec_rad"])
            Rawx.extend(uvh5debug["matrix_h"])

            # Update non-scalar Header datasets
            uvh5hdr["ant_1_array"][(1:nbls).+nblts] = a1s
            uvh5hdr["ant_2_array"][(1:nbls).+nblts] = a2s
            uvh5hdr["uvw_array"][:,(1:nbls).+nblts] = bl_uvws
            uvh5hdr["time_array"][(1:nbls).+nblts] = timejd
            uvh5hdr["lst_array"][(1:nbls).+nblts] = lst
            uvh5hdr["integration_time"][(1:nbls).+nblts] = tau

            uvh5debug["dut1"][1+ntimes] = ut1utcsec
            uvh5debug["ra_rad"][1+ntimes] = ra_rad
            uvh5debug["dec_rad"][1+ntimes] = dec_rad
            uvh5debug["matrix_h"][:, :, 1+ntimes] = matrix_h

            # Update Data datasets
            # Flag data if all subints were zero
            uvh5data["flags"][:, :, (1:nbls).+nblts] = subints_per_dump == subints_zero ? 1 : 0
            # Set nsamples to fraction on non-zero subints
            uvh5data["nsamples"][:, :, (1:nbls).+nblts] = (subints_per_dump - subints_zero) / subints_per_dump
            subints_zero = 0

            idx = 1
            for (xgpuidx, blidx, polidx, isauto, needsconj) in inpair_maps
                idx += 1

                visdata[polidx, :, blidx] = (
                    needsconj ? conj.(matrix_h[xgpuidx, 1:nfreqs])
                              :       matrix_h[xgpuidx, 1:nfreqs])
                # If cross-pol autocorrelation
                if isauto && polidx == 2
                  # Output redundant cross-pol autocorrelation
                  # We could cheat and output conjugate of current xgpuidx
                  # spectrum, but instead we use some inside knowledge that
                  # should be publicized in XGPU documentation that the
                  # redundant cross-pol is at xgpuidx-1.
                  visdata[polidx+1, :, blidx] = (
                      needsconj ? conj.(matrix_h[xgpuidx-1, 1:nfreqs])
                                :       matrix_h[xgpuidx-1, 1:nfreqs])
                end
            end
            uvh5data["visdata"][:, :, (1:nbls).+nblts] = visdata

            # Update scalar datasets in Header
            nblts += nbls
            ntimes += 1
            write(uvh5hdr["Nblts"], nblts)
            write(uvh5hdr["Ntimes"], ntimes)
          end # if dump time
        end # for subint in block

        # Increment pktidx
        pktidx += piperblk
      end # block loop
      # Break out of "while read loop" after first dump
      break
    end # while read header
  end # do open raw

  @info "output $ntimes integrations"
  if subints_per_dump > 0
    @info "$(subints_sent%subints_per_dump)/$subints_per_dump leftover"
  end

  # Dump debug datasets to test headers
  open("matrix_h.h", "w") do matfio
    matrix_h = matrix_h[:,:,1] # take first Ntimes dim
    write(matfio, "#ifndef _MATRIX_H_H\n")
    write(matfio, "#define _MATRIX_H_H\n\n")
    write(matfio, "#define MATRIX_H_RANK ", string(length(size(matrix_h))+1), "\n")
    write(matfio, "// MATRIX_H_DIMS [Nfreqs, xGPU_prods_out, real_imag]\n")
    write(matfio, "const int MATRIX_H_DIMS[MATRIX_H_RANK] = {", join(reverse(size(matrix_h)), ","), ",2};\n\n")
    
    write(matfio, "const int MATRIX_H[] = {\n\t")
    map(i->write(matfio, replace(replace(string(i), r"\s\+?([-]?)\s*"=>s", \1"), "im"=>",\n\t")), matrix_h)
    write(matfio, "};\n\n#endif // _MATRIX_H_H\n")
  end
  open("visdata.h", "w") do vdfio
    write(vdfio, "#ifndef _VISDATA_H\n")
    write(vdfio, "#define _VISDATA_H\n\n")
    write(vdfio, "#define VISDATA_RANK ", string(length(size(visdata))+1), "\n")
    write(vdfio, "// VISDATA_DIMS [Nblts, Nfreqs, Npols, real_imag]\n")
    write(vdfio, "const int VISDATA_DIMS[VISDATA_RANK] = {", join(reverse(size(visdata)), ","), ",2};\n\n")
    
    write(vdfio, "const int VISDATA[] = {\n\t")
    map(i->write(vdfio, replace(replace(string(i), r"\s\+?([-]?)\s*"=>s", \1"), "im"=>",\n\t")), visdata)
    write(vdfio, "};\n\n#endif // _VISDATA_H\n")
  end
  
  open("uvw_phased.h", "w") do uvwfio
    ut1utcsec = uvh5debug["dut1"][ntimes]
    ra_rad = uvh5debug["ra_rad"][ntimes]
    dec_rad = uvh5debug["dec_rad"][ntimes]
    write(uvwfio, "#ifndef _UVW_PHASED_H\n")
    write(uvwfio, "#define _UVW_PHASED_H\n\n")
    write(uvwfio, "const double hour_angle_rad =  ", string(hob), ";\n")
    write(uvwfio, "const double declination_rad =  ", string(dob), ";\n")
    write(uvwfio, "const double timemjd =  ", string(timejd), ";\n")
    write(uvwfio, "const double dut1 =  ", string(ut1utcsec), ";\n")
    write(uvwfio, "const double ra_rad =  ", string(ra_rad), ";\n")
    write(uvwfio, "const double dec_rad =  ", string(dec_rad), ";\n\n")

    write(uvwfio, "#define ENU_RANK ", string(length(size(antpos_enu))), "\n")
    write(uvwfio, "// ENU_DIMS [Nants_data, 3]\n")
    write(uvwfio, "const double ENU_DIMS[ENU_RANK] = {", join(reverse(size(antpos_enu)), ","), "};\n")
    write(uvwfio, "const double ENU[] = {\n\t")
    map(enu->write(uvwfio, join(map(string, enu), ", ")*",\n\t"), antpos_enu[1:3, i] for i in 1:size(antpos_enu, 2))
    write(uvwfio, "};\n\n")

    # write(uvwfio, "#define BL_RANK 2\n")
    # write(uvwfio, "// BL_DIMS [Nbls 2]\n")
    # write(uvwfio, "const int BL_DIMS[BL_RANK] = {", string(length(bl_array)), ",2};\n")
    # write(uvwfio, "const int BL[] = {\n\t")
    # map(bl->write(uvwfio, join(map(string, bl.-1), ", ")*",\n\t"), bl_array)
    # write(uvwfio, "};\n\n")

    write(uvwfio, "#define ANTENNA_ARRAY_RANK 1\n")
    write(uvwfio, "// ANTENNA_ARRAY_DIMS [Nbls]\n")
    write(uvwfio, "const int ANTENNA_ARRAY_DIMS[ANTENNA_ARRAY_RANK] = {", string(length(a1s)), "};\n")
    write(uvwfio, "const int ANTENNA_1_ARRAY[] = {\n\t")
    map(a1->write(uvwfio, string(a1)*",\n\t"), a1s)
    write(uvwfio, "};\n\n")
    write(uvwfio, "const int ANTENNA_2_ARRAY[] = {\n\t")
    map(a2->write(uvwfio, string(a2)*",\n\t"), a2s)
    write(uvwfio, "};\n\n")

    write(uvwfio, "#define UVW_RANK ", string(length(size(bl_uvws))), "\n")
    write(uvwfio, "// UVW_DIMS [Nblts, 3]\n")
    write(uvwfio, "const double UVW_DIMS[UVW_RANK] = {", join(reverse(size(bl_uvws)), ","), "};\n")
    write(uvwfio, "const double UVW[] = {\n\t")
    map(uvw->write(uvwfio, join(map(string, uvw), ", ")*",\n\t"), bl_uvws[1:3, i] for i in 1:size(bl_uvws, 2))
    write(uvwfio, "};\n\n#endif // _UVW_PHASED_H\n")
  end

  # Close uvh5
  if uvh5 !== nothing
    close(uvh5debug)
    close(uvh5)
  end

  # Free XGPU resources
  if xgpuctx !== nothing
    xgpuFree(xgpuctx)
  end

end # uvh5debug_fromraw()



telinfo = "/home/sonata/telinfo_ata.yml"
obsinfo = "/home/sonata/scripts/obsinfo.yml"

inttime_per_dump = 0.016384

rawstem = "/home/sonata/dev/hdf5_tests/guppi_59596_67211_481811_3c380_0001.2block"
uvh5file = "/home/sonata/dev/hdf5_tests/guppi_59596_67211_481811_3c380_0001.debug.uvh5"
cd(@__DIR__)
if isfile(uvh5file)
  @info "removing existing output file $uvh5file"
  rm(uvh5file)
end

# Uses Julia >=1.5 bare keyword argument feature!
uvh5debug_fromraw(rawstem, uvh5file; telinfo, obsinfo, inttime_per_dump)