using Images
using CUDA

function loadholo(path::String)
    out = Float32.(channelview(Gray.(load(path))))
end

function CuTransSqr!(Plane::CuArray, datLen::Int, wavLen::Float32, dx::Float32)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        Plane[x,y] = 1.0 - ((x-datLen/2)*wavLen/datLen/dx)^2 - ((y-datLen/2)*wavLen/datLen/dx)^2
    end
end

function CuTransFunc!(Plane::CuArray, d_sqrPart::CuArray, z0::Float32, wavLen::Float32, datLen::Int)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        Plane[x,y] = exp(2im*pi*(z0)/wavLen*sqrt(d_sqrPart[x,y]))
    end
end

function CuUpdateImposed!(imposed::CuArray, input::CuArray, datLen::Int)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y
    if x <= datLen && y <= datLen
        if input[x,y] < imposed[x,y]
            imposed[x,y] = input[x,y]
        end
    end
end


function phaseRetrieval!(Plane::CuArray, img1::CuArray, img2::CuArray, m::Int, dz::Float32,)
    
end