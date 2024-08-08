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


# Load the original image
p1 = Gray.(load("Mri_Brain_Flair_Axial.jpg"))
p2 = zeros(Float64,size(p1,1),size(p1,2))
p3 = zeros(Float64,size(p1,1),size(p1,2))
p2 = load_image(p1, p2)

# Calling Bspline interpolation function
itp = interpolate(p2, BSpline(Cubic(Reflect(OnGrid()))))

# Desired translation/rotaion parameters
x = [-15.0, 3.0, 0.3]

N = size(p1,1)
M = size(p1,2)

# Manipulate each pixel value in the image
for i in 1:size(p1, 1)
    for j in 1:size(p1, 2)
        # Scale the pixel value
        p3[i, j] = itp(min(max(cos(x[3])*(i-N/2) + sin(x[3])*(j-M/2) + x[1] + N/2, 1), size(p1, 1)), 
        min(max(-sin(x[3])*(i - N/2) + cos(x[3])*(j-M/2) + x[2] + M/2, 1), size(p1, 2)))
                
    end
end

# Save the manipulated image with your dsired name and path
save("Mri_Brain_Flair_Axial2.jpg", p3)
