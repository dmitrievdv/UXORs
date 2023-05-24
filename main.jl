using FITSIO
using Plots
using Printf
using DelimitedFiles
using Dates

function findbkg(cut, aperture_size)
    cut_width = size(cut)[1]
    cut_height = size(cut)[2]
    aperture_flux = zeros(cut_width-2aperture_size, cut_height-2aperture_size)
    for i=1+aperture_size:cut_width-aperture_size
        width_aperture = i-aperture_size:i+aperture_size
        for j=1+aperture_size:cut_height-aperture_size
            height_aperture = j-aperture_size:j+aperture_size
            aperture_flux[i-1,j-1] = sum(cut[width_aperture, height_aperture])
        end
    end
    min_index = findmin(aperture_flux)[2] + CartesianIndex(1,1)
    return min_index[1]-aperture_size:min_index[1]+aperture_size, min_index[2]-aperture_size:min_index[2]+aperture_size
end

negativetonan(x) = x < 0 ? NaN : x

function findstar(cut, aperture_size)
    cut_width = size(cut)[1]
    cut_height = size(cut)[2]
    aperture_flux = zeros(cut_width-2aperture_size, cut_height-2aperture_size)
    for i=1+aperture_size:cut_width-aperture_size
        width_aperture = i-aperture_size:i+aperture_size
        for j=1+aperture_size:cut_height-aperture_size
            height_aperture = j-aperture_size:j+aperture_size
            aperture_flux[i-1,j-1] = sum(cut[width_aperture, height_aperture])
        end
    end
    max_index = findmax(aperture_flux)[2] + CartesianIndex(1,1)
    return max_index[1]-aperture_size:max_index[1]+aperture_size, max_index[2]-aperture_size:max_index[2]+aperture_size
end

function getlc(star, tess_cut_fits)
    sector = parse(Int, tess_cut_fits[7:10])
    fits = FITS("$star/$tess_cut_fits")
    cuts = negativetonan.(read(fits[2], "FLUX"))
    cuts_shape = size(cuts)
    cut_width = cuts_shape[1]
    cut_height = cuts_shape[2]
    N_cuts = cuts_shape[3]
    N_bkg = N_cuts ÷ 4

    jd = read(fits[2], "TIME") .+ 2457000

    plt = heatmap(cuts[:,:,N_bkg]', aspect_ratio = :equal)
    

    aperture_size = 1
    center_px_x = (cut_width+1) ÷ 2
    center_px_y = (cut_height+1) ÷ 2
    aperture = (center_px_x-2aperture_size:center_px_x+2aperture_size, center_px_y-2aperture_size:center_px_y+2aperture_size)
    aperture = findstar(cuts[aperture..., N_bkg], aperture_size)
    aperture = aperture[1] .+ (center_px_x - 2aperture_size - 1), aperture[2] .+ (center_px_y - 2aperture_size - 1)
    # aperture = (center_px_x-aperture_size:center_px_x+aperture_size, center_px_y-aperture_size:center_px_y+aperture_size)
    bkg_aperture = findbkg(cuts[:,:,N_bkg], aperture_size)
    raw_flux = reshape(sum(cuts[aperture..., :], dims = (1,2)), N_cuts)
    bkg_flux = reshape(sum(cuts[bkg_aperture..., :], dims = (1,2)), N_cuts)
    flux = negativetonan.(raw_flux .- bkg_flux)

    bkg_aperture_rect_x = [bkg_aperture[1][1] - 0.5, bkg_aperture[1][1] - 0.5, 
                           bkg_aperture[1][end] + 0.5, bkg_aperture[1][end] + 0.5, 
                           bkg_aperture[1][1] - 0.5]

    bkg_aperture_rect_y = [bkg_aperture[2][1] - 0.5, bkg_aperture[2][end] + 0.5, 
                           bkg_aperture[2][end] + 0.5, bkg_aperture[2][1] - 0.5, 
                           bkg_aperture[2][1] - 0.5]

    plot!(plt, bkg_aperture_rect_x, bkg_aperture_rect_y, lw = 2, label = "bkg", rightmargin = 20px)


    aperture_rect_x = [aperture[1][1] - 0.5, aperture[1][1] - 0.5, 
                           aperture[1][end] + 0.5, aperture[1][end] + 0.5, 
                           aperture[1][1] - 0.5]

    aperture_rect_y = [aperture[2][1] - 0.5, aperture[2][end] + 0.5, 
                           aperture[2][end] + 0.5, aperture[2][1] - 0.5, 
                           aperture[2][1] - 0.5]

    plot!(plt, aperture_rect_x, aperture_rect_y, lw = 2, label = "star", rightmargin = 20px)

    savefig(plt, "$star/$sector-map.pdf")
    savefig(plt, "maps/$star-$sector-map.pdf")

    open("$star/$sector-rawlc.dat", "w") do io
        println(io, "#jd flux raw_flux bkg_flux")
        for n = 1:N_cuts
            @printf(io, "%15.6f %12.3e %12.3e %12.3e\n", jd[n], flux[n], raw_flux[n], bkg_flux[n])
        end
    end

    return Float64.(jd), Float64.(flux)
    # uniform_jd_step = round(minimum(jd[2:end] - jd[1:end-1])*24*60)/24/60
    # jd_start = jd[1]
    # jd_end = jd[end]
end

function uniformlc(jd, flux, ε_jd = 1.5, ε_flux = 0.1)
    uniform_jd = Float64[]
    uniform_flux = Float64[]
    N_samples = length(jd)
    jd_step = jd[2] - jd[1]
    flux_step = flux[2] - flux[1]
    for i = 2:N_samples
        new_jd_step = jd[i] - jd[i-1]
        new_flux_step = flux[i] - flux[i-1]
        if (new_jd_step/jd_step > ε_jd) | (jd_step/new_jd_step > ε_jd)
            push!(uniform_jd, jd[i])
            push!(uniform_flux, NaN)
        elseif (abs(new_flux_step)/flux[i] > ε_flux)
            push!(uniform_jd, jd[i])
            push!(uniform_flux, NaN)
            if abs(flux_step + new_flux_step)/flux[i] < ε_flux
                push!(uniform_jd, jd[i])
                push!(uniform_flux, flux[i])
            end
        else
            push!(uniform_jd, jd[i])
            push!(uniform_flux, flux[i])
        end
        jd_step = new_jd_step
        flux_step = new_flux_step
    end
    
    return uniform_jd, uniform_flux
end

stars = string.(readdlm("UXORs.data")[:, 1])

for star in stars
    println(star)
    for fits_file in readdir(star)
        if fits_file[1:4] != "tess"
            continue
        end
        println(fits_file)
        # fits_file="tess-s0032-1-2_76.123670_-3.787760_15x15_astrocut.fits"
        sector = parse(Int, fits_file[7:10])

        jd, flux = getlc(star, fits_file)
        jd, flux = uniformlc(jd, flux)
        
        mag = -2.5*log10.(flux) .+ 20.44

        open("$star/$sector-lc.dat", "w") do io
            println(io, "#jd flux")
            for n = 1:length(jd)
                @printf(io, "%15.6f %12.3e %12.3f\n", jd[n], flux[n], mag[n])
            end
        end

        int_dates = collect(ceil(Int, jd[1]/5)*5:5:floor(Int, jd[end]/5)*5)
        string_dates = Dates.format.(julian2datetime.(int_dates), "d u Y")
        plt = plot(jd, mag, xticks = (int_dates, string_dates), label = false, rightmargin = 15px, yflip = true)
        savefig(plt, "plots/$star-$sector.pdf")
        savefig(plt, "$star/$sector-lc.pdf")
    end
end
# plot(jd, flux)
# uniformlc(jd, flux)


