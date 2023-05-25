using DelimitedFiles
using Plots
using Plots.PlotMeasures
using FFTW

function splitlc(jd, flux)
    split_jd = Vector{Float64}[]
    split_flux = Vector{Float64}[]
    block_start = 1
    block_end = 1
    for i = eachindex(jd)
        if isnan(flux[i])
            block_end = i-1
            if (block_end - block_start > 100)
                push!(split_jd, jd[block_start:block_end])
                push!(split_flux, flux[block_start:block_end])
            end
            block_start = i+1
        end
    end
    push!(split_jd, jd[block_start:end])
    push!(split_flux, flux[block_start:end])
    return split_jd, split_flux
end

function freq(jd)
    N_samples = length(jd)
    sampling_rate = 1/(jd[2]- jd[1])
    freq = fftfreq(N_samples, sampling_rate)
    return freq
end

function finddisp(freq, abs_fft)
    freq_step = freq[2] - freq[1]
    norm = sum(abs_fft[2:end])*freq_step
    return sum(freq[2:end] .^ 2 .* abs_fft[2:end] /norm)*freq_step
end

star = "COOri"
lc_file = "$star/lc.dat"
lc_data = Float64.(readdlm(lc_file)[2:end,:])

jd = lc_data[:,1]
flux = lc_data[:,2]

split_jd, split_flux = splitlc(jd, flux)
# split_flux[2] = split_flux[2] .- 1000
freqs = freq.(split_jd)
ffts = fft.(split_flux)
abs_ffts = [abs.(fft) for fft in ffts]

finddisp.(freqs, abs_ffts)

plt_fft = plot(fftshift.(freqs), fftshift.(abs_ffts), xlims = (-5,5))
# plt_lc = plot(split_jd, split_flux)