using DelimitedFiles
using Plots
using Plots.PlotMeasures
using FFTW
using NFFT
using LinearRegression
using LombScargle

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
    N = length(jd)
    if (N - block_start > 100)
        push!(split_jd, jd[block_start:end])
        push!(split_flux, flux[block_start:end])
    end
    return split_jd, split_flux
end

function freq(jd)
    N_samples = length(jd)
    sampling_rate = 1/(jd[2]- jd[1])
    freq = fftfreq(N_samples, sampling_rate)
    return freq
end

function findfftdisp(freq, fft_dist)
    freq_step = freq[2] - freq[1]
    return sum(freq .^ 2 .* fft_dist)*freq_step
end

function absft(jd, fft)
    sampling_step = jd[2] - jd[1]
    return abs.(fft)*sampling_step
end

function fftdist(freq, abs_fft)
    freq_step = freq[2] - freq[1]
    norm = sum(abs_fft[1:end])*freq_step
    return abs_fft/norm
end

function detrendlinearregresion(jd, flux)
    jd_start = jd[1]; jd_end = jd[end]
    times = jd .- jd_start
    lr = linregress(times, flux)
    a, b = coef(lr)
    detrended_flux = flux - a*times .- b
end

function detrend(jd, flux)
    N_samples = length(jd)
    jd_start = jd[1]; jd_end = jd[end]
    jd_range = jd_end - jd_start
    flux_start = flux[1]; flux_end = flux[end]
    flux_change = flux_end - flux_start
    detrended_flux = flux .- (flux_start .+ flux_change/jd_range * (jd .- jd_start))
    return detrended_flux
end

function detrendandextend(jd, flux)
    N_samples = length(jd)
    jd_start = jd[1]; jd_end = jd[end]
    jd_range = jd_end - jd_start
    flux_start = flux[1]; flux_end = flux[end]
    flux_change = flux_end - flux_start
    detrended_flux = flux .- (flux_start .+ flux_change/jd_range * (jd .- jd_start))
    extended_flux = vcat(detrended_flux, -reverse(detrended_flux)[2:end-1])
    # extended_jd = vcat(jd .- jd_start, jd_end .+ (jd[2:end-1] .- jd_start))
    return extended_flux
end

function extendjd(jd)
    jd_start = jd[1]; jd_end = jd[end]
    jd_range = jd_end - jd_start
    extended_jd = vcat(jd .- jd_start, jd_range .+ (jd[2:end-1] .- jd_start))
    return extended_jd
end

function derivfft(freq, fft)
    return 1im*2π*freq .* fft
end

function trendderivative(jd, flux)
    jd_start = jd[1]; jd_end = jd[end]
    jd_range = jd_end - jd_start
    flux_start = flux[1]; flux_end = flux[end]
    flux_change = flux_end - flux_start
    return flux_change/jd_range
end

function deextend(array)
    N_samples = length(array) ÷ 2 + 1
    return array[1:N_samples]
end

function pltderiv(n; kwargs...)
    plt = plot(split_jd[n], split_flux[n] .- split_flux[n][1]; kwargs...)
    # plot!(plt, split_jd[n][1:end-1], (split_flux[n][2:end] - split_flux[n][1:end-1]) ./ (split_jd[n][2:end] - split_jd[n][1:end-1]))
    plot!(split_jd[n], split_deriv[n])
end

function derivdisp(deriv)
    √(sum(deriv .^ 2)/length(deriv))
end

smoothfft(freq, fft) = [abs(freq[i]) > 8 ? 1e-50im : fft[i] for i = 1:length(freq)]

samplingstep(jd) = jd[2]-jd[1]

function getffts(split_jd, split_flux)
    extended_fluxes = detrendlinearregresion.(split_jd, split_flux)
    extended_jds = split_jd
    # plt_lc = plot(extended_jds, extended_fluxes)
    # split_flux[2] = split_flux[2] .- 1000
    freqs = freq.(extended_jds)
    ffts = fft.(extended_fluxes)
    return freqs, ffts
end

function getffts(star)
    lc_file = "$star/lc.dat"
    lc_data = try
        Float64.(readdlm(lc_file)[2:end,:])
    catch
        return -1
    end
    jd = lc_data[:,1]
    flux = lc_data[:,3]

    split_jd, split_flux = splitlc(jd, flux)
    extended_fluxes = detrend.(split_jd, split_flux)
    extended_jds = split_jd
    # plt_lc = plot(extended_jds, extended_fluxes)
    # split_flux[2] = split_flux[2] .- 1000
    freqs = freq.(extended_jds)
    ffts = fft.(extended_fluxes)
    return freqs, ffts
end

function splitlc(star)
    lc_file = "$star/lc.dat"
    lc_data = try
        Float64.(readdlm(lc_file)[2:end,:])
    catch
        return -1
    end
    jd = lc_data[:,1]
    flux = lc_data[:,3]

    splitlc(jd, flux)
end

function getderivdisp(star)
# star = "BMAnd"
    lc_file = "$star/lc.dat"
    lc_data = try
        Float64.(readdlm(lc_file)[2:end,:])
    catch
        return -1
    end
    jd = lc_data[:,1]
    flux = lc_data[:,3]

    split_jd, split_flux = splitlc(jd, flux)
    extended_fluxes = detrendandextend.(split_jd, split_flux)
    extended_jds = extendjd.(split_jd)
    plt_lc = plot(extended_jds, extended_fluxes)
    # split_flux[2] = split_flux[2] .- 1000
    freqs = freq.(extended_jds)
    ffts = fft.(extended_fluxes)
    smooth_ffts = smoothfft.(freqs, ffts)
    deriv_ffts = derivfft.(freqs, smooth_ffts)
    deriv_extended_fluxes = [real.(ifft(fft)) for fft in deriv_ffts]
    split_deriv = [deextend(deriv_extended_fluxes[i]) .+ trendderivative(split_jd[i], split_flux[i]) for i = 1:length(split_jd)]
    deriv_jd = deextend.(extended_jds)
    abs_fts = absft.(extended_jds, smooth_ffts)
    fft_dists = fftdist.(freqs, abs_fts)

    disp = derivdisp.(split_deriv)
    max = [maximum(abs.(deriv)) for deriv in split_deriv]
    println(disp)
    
    plt_fft = plot(fftshift.(freqs), fftshift.(fft_dists), xlims = (-5,5), yrange = (-0.5,3), title = "$star")
    plt_deriv = plot(deriv_jd, split_deriv)
    return sum(disp)/length(disp)
end

function joincloselc(jds, lcs; Δt = 2.0)
    joined_jds = Vector{Float64}[]
    joined_lcs = Vector{Float64}[]
    push!(joined_jds, deepcopy(jds[1]))
    push!(joined_lcs, deepcopy(lcs[1]))
    for (jd, lc) in zip(jds[2:end], lcs[2:end])
        if jd[1] - joined_jds[end][end] < Δt
            append!(joined_jds[end], deepcopy(jd))
            append!(joined_lcs[end], deepcopy(lc))
        else
            push!(joined_jds, deepcopy(jd))
            push!(joined_lcs, deepcopy(lc))
        end
    end
    return joined_jds, joined_lcs
end

function getNFFT(jd, lc; N2 = 1000)
    N = 2*N2
    k = (jd .- jd[1])
    t_end = k[end]
    k = k / 10t_end .- 0.5
    J = length(k)
    f = deepcopy(lc)
    p = plan_nfft(k, N, reltol = 1e-9, σ = 2.0)
    freqs = [-N2:N2-1;] / 10t_end
    return freqs, adjoint(p) * f
end


star = "WWVul"


jds, lcs = splitlc(star)
detrended_lcs = detrendlinearregresion.(jds, lcs) 

freqs, ffts = getffts(jds, lcs)
joined_jds, joined_lcs = joincloselc(jds, lcs)
joined_jds, joined_detrended_lcs = joincloselc(jds, detrended_lcs, Δt = 2)

nfft_freqs = Vector{Float64}[]
nffts = Vector{ComplexF64}[]
for (jd, lc) in zip(joined_jds, joined_detrended_lcs)
    nfft_freq, nfft = getNFFT(jd, lc)
    push!(nfft_freqs, nfft_freq)
    push!(nffts, nfft)
end

nfft_phases = zeros(length(nffts))
nfft_periods = zeros(length(nffts))
for i = eachindex(nfft_phases)
    i_max = findmax(abs.(nffts[i]) .^ 2)[2]
    nfft_periods[i] = 1 / abs(nfft_freqs[i][i_max])
    
    nfft_phases[i] = angle(nffts[i][i_max]) 
    # if nfft_phases[i] < 0.0 
    #     nfft_phases[i] += π
    # end
end

# joined_detrended_lcs = detrendlinearregresion.(joined_jds, joined_lcs)
periodograms = lombscargle.(joined_jds, joined_detrended_lcs)

inds = [1,2]


plt_joined = plot(legend = false, yflip = true)
for (jd, lc) in zip(joined_jds, joined_lcs)
    plot!(plt_joined, jd .- jd[1], lc)
    # scatter!(plt_joined, jd .- jd[1], lc, ms = 1, mc = :black)
end
plt_joined

plt_joined_detrended = plot(legend = false, yflip = true)
for (jd, lc) in zip(joined_jds, joined_detrended_lcs)
    plot!(plt_joined_detrended, jd .- jd[1], lc)
    # scatter!(plt_joined_detrended, jd .- jd[1], lc, ms = 1, mc = :black)
end
plt_joined_detrended

# nfft_period = sum(nfft_periods)/length(nfft_periods)
plt_joined_phased= plot(legend = false, yflip = true, xlims = (0,25))
for (jd, lc, ph, per) in zip(joined_jds, joined_detrended_lcs, nfft_phases, nfft_periods)
    plot!(plt_joined_phased, jd .- jd[1] .+ ph/2π*per, lc)
    # scatter!(plt_joined_detrended, jd .- jd[1], lc, ms = 1, mc = :black)
end
plt_joined_phased

plt_periodograms = plot(legend = false, xlims = (0,1))
for periodogram in periodograms
    plot!(plt_periodograms, freqpower(periodogram)...)
end
plt_periodograms