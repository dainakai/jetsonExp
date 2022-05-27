using CUDA

function elementwiseProc(array1,array2,datLen)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y

    if x <= datLen && y <= datLen
        array1[y,x] = array1[y,x]^array2[y,x]
    end
    return nothing
end

function CuGetCrossCor!(corArray::CuDeviceArray{Float32,2},img1::CuDeviceArray{Float32,2},img2::CuDeviceArray{Float32,2},gridNum::Int64,srchSize::Int64,intrSize::Int64,gridSize::Int64)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y


    if x <= (srchSize-intrSize+1)*(gridNum-1) && y <= (srchSize-intrSize+1)*(gridNum-1)
        gridIdxx = div(x-1,(srchSize-intrSize+1))+1
        gridIdxy = div(y-1,(srchSize-intrSize+1))+1
        idxx  = x - (gridIdxx-1)*(srchSize-intrSize+1)
        idxy  = y - (gridIdxy-1)*(srchSize-intrSize+1)

        a1::Int64 = gridIdxy*gridSize-div(intrSize,2)
        a2::Int64 = gridIdxx*gridSize-div(intrSize,2)
        b1::Int64 = (gridIdxy-1)*gridSize+idxy-1
        b2::Int64 = (gridIdxx-1)*gridSize+idxx-1

        meanA::Float32 = 0.0
        meanB::Float32 = 0.0
        num::Float32 = 0.0
        denomA::Float32 = 0.0
        denomB::Float32 = 0.0

        for i in 1:intrSize
            for j in 1:intrSize
                meanA += img1[a1+i,a2+j]
                meanB += img2[b1+i,b2+j]
            end
        end
        meanA /= Float32(intrSize^2)
        meanB /= Float32(intrSize^2)

        for i in 1:intrSize
            for j in 1:intrSize
                num += (img1[a1+i,a2+j]-meanA)*(img2[b1+i,b2+j]-meanB)
                denomA += (img1[a1+i,a2+j]-meanA)^2
                denomB += (img2[b1+i,b2+j]-meanB)^2
            end
        end

        corArray[y,x] = num/(CUDA.sqrt(denomA)*CUDA.sqrt(denomB))
    end
    return nothing
end

function main()
    datLen = 1024*2^2

    array1 = CUDA.rand(datLen,datLen)
    array2 = CUDA.rand(datLen,datLen)

    blockSize = 16
    threads1 = (blockSize,blockSize)
    blocks1 = (cld(datLen,blockSize),cld(datLen,blockSize))

    @cuda threads = threads1 blocks = blocks1 elementwiseProc(array1,array2,datLen)
    
    hostA1 = Array(array1)
    display(hostA1)
    return nothing
end

main()