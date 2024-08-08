using Optimization, ForwardDiff, ReverseDiff, Zygote, OptimizationOptimJL, Plots
using Images, ImageIO, Interpolations

# Loading gray-scale images and extracting their pixel values
function load_image(p1,p2)
    for r in 1:size(p1, 1)
        for c in 1:size(p1, 2)
            p2[r,c] = Float64(p1[r,c].val)
        end
    end
    return p2
end


# This function applies the obtained parameters at each level of the pyramid to make a correct image
function transfer(p11, t1, t2, t3, i, output_dir)

    itp = interpolate(p11, BSpline(Cubic(Reflect(OnGrid()))))

    x = [t1, t2, t3]

    N = size(p11,1)
    M = size(p11,2)
    # Scale each pixel value in the image
    for i in 1:size(p11, 1)
        for j in 1:size(p11, 2)
            # Scale the pixel value
            p11[i, j] = itp(min(max(cos(x[3])*(i-N/2) + sin(x[3])*(j-M/2) + x[1] + N/2, 1), size(p11, 1)), 
            min(max(-sin(x[3])*(i - N/2) + cos(x[3])*(j-M/2) + x[2] + M/2, 1), size(p11, 2)))

        end
    end

    #save transformed images
    save("$(output_dir)transformed$(i)1.png", p11)

end


# This is the main core of the code - Optimization function 
function optimize(p3, p4, t1, t2, t3)
    itp = interpolate(p3, BSpline(Cubic(Reflect(OnGrid()))))

    x0 = [t1, t2, t3]
    p = [itp, p4]

    # Defining the objective function
    function myfunc(x, p)
        itp, p4 = p
        N, M = size(p4)
        res = 0.0
        cos_x3 = cos(x[3])
        sin_x3 = sin(x[3])
        half_N = N / 2
        half_M = M / 2

        for i in 1:N
            for j in 1:M
                xi = cos_x3 * (i - half_N) + sin_x3 * (j - half_M) + x[1] + half_N
                yi = -sin_x3 * (i - half_N) + cos_x3 * (j - half_M) + x[2] + half_M
                xi = clamp(xi, 1, N)
                yi = clamp(yi, 1, M)

                value_p1 = itp(xi, yi)
                value_p2 = p4[i, j]
                diff = (value_p1 - value_p2)^2
                res += diff
            end
        end
        return res
    end

    # Create the optimization function and problem
    f = OptimizationFunction(myfunc, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(f, x0, p)
    sol = solve(prob, BFGS(), maxiters=10000, show_trace=false)
    t1, t2, t3 = sol.u

    return t1, t2, t3
end

# Levels of the Gaussian Pyramid
n_scales = 5

# Ratio between the size of the images in two consecutive pyramid levels
downsample = 2

# Size of the Gaussian kernel
sigma = 2.0


# Obtaining the Gaussian pyramid for the reference and transformed images
pyramid1 = gaussian_pyramid(Gray.(load("Mri_Brain_Flair_Axial2.jpg")), n_scales, downsample, sigma)
pyramid2 = gaussian_pyramid(Gray.(load("Mri_Brain_Flair_Axial.jpg")), n_scales, downsample, sigma)


# Path to save your generated images
output_dir = "C:/Users/user/Desktop/julia_codes/pyramid/"
num_levels = 5 # normally, should be equal to n_scales


# Loop to save images and combine them
for i in 1:num_levels
    save("$(output_dir)downsampled$(i)1.png", pyramid1[i])
    save("$(output_dir)downsampled_reference$(i)2.png", pyramid2[i])
    
    combined_image = hcat(pyramid1[i], pyramid2[i])
    save("$(output_dir)combined$(i).png", combined_image)
end


# Initializing zero matrices
for i in 1:num_levels
    eval(Meta.parse("global p$(i)1 = zeros(Float64, size(pyramid1[$i], 1), size(pyramid1[$i], 2))"))
    eval(Meta.parse("global p$(i)2 = zeros(Float64, size(pyramid2[$i], 1), size(pyramid2[$i], 2))"))
end


# Loading the images into the matrices
for i in 1:num_levels
    eval(Meta.parse("global p$(i)1 = load_image(pyramid1[$i], p$(i)1)"))
    eval(Meta.parse("global p$(i)2 = load_image(pyramid2[$i], p$(i)2)"))
end


# Intialization of the translation/rotaion parameters
t1, t2, t3 = 0.0, 0.0, 0.0


# Optimizer loop through each level
for i in 1:num_levels
    eval(Meta.parse("global t1, t2, t3 = optimize(p$(num_levels+1-i)1, p$(num_levels+1-i)2, t1, t2, t3)"))
    eval(Meta.parse("transfer(p$(num_levels+1-i)1, t1, t2, t3, $(num_levels+1-i), output_dir)"))
    println("Shift in i direction: ", t1, " Shift in j direction: ", t2, " angle: ", t3, '\n')
    global t1 = t1 * downsample
    global t2 = t2 * downsample
end





#= Uncomment in case of timing inspections but you shoud use "BenchmarkTools" package and call
@btime begin
end
In the beggining and the end of the part you want to do timing inspections

print(" time at level 5: ", @time optimize(p51, p52, t1, t2, t3), " time at level 4: ", @time optimize(p41, p42, t1, t2, t3), 
 " time at level 3: ", @time optimize(p31, p32, t1, t2, t3), " time at level 2: ", @time optimize(p21, p22, t1, t2, t3),
 " time at level 1: ", @time optimize(p11, p12, t1, t2, t3))
 =#