"""
2枚の画像をPIVしてベクトル場と一緒にGLMakieで表示
"""

# julia -t 4 とか。
using Images
# using GLMakie
using CairoMakie
using StatsBase
using CUDA
# include("functions.jl")

"""
修正コレスキー分解

対象行列Aを A = LDL^T に変換します。

戻り値行列の対角成分は D の逆数、左下成分は L に一致します。
"""
function choleskyDecom(Matrix)
    ndims(Matrix) == 2 ? 1 : error("choleskyDecom : Input Matrix is not in 2 dims.")
    n = size(Matrix)[1]
    n == size(Matrix)[2] ? 1 : error("choleskyDecom : Input Matrix is not square.")
    # println("Input : $(n)*$(n) matrix")
    matA = copy(Matrix)
    vecw = Array{Float64}(undef,n)
    for j in 1:n
        for i in 1:j-1
            vecw[i] = matA[j,i]
            for k in 1:i-1
                vecw[i] -= matA[i,k] * vecw[k]
            end
            matA[j,i] = vecw[i] * matA[i,i]
        end
        t = matA[j,j]
        for k in 1:j-1
            t -= matA[j,k] * vecw[k]
        end
        matA[j,j] = 1.0/t
    end
    # matA の対角成分は matD の逆数。matA_ji (i<j) は matL_ji。
    return matA
end

"""
連立方程式ソルバー

修正コレスキー分解で取得した行列 `chlskyMat` およびヤコビアン `yacob` 、エラーベクトル `errorArray`から

`chlskyMat` x = - `yacob` `errorArray`

を解き、戻り値として返します。
"""
function simEqSolver(chlskyMat,yacob,errorArray)
    vecB = -transpose(yacob)*errorArray
    n = size(vecB)[1]

    matL = chlskyMat
    

    vecX = zeros(n)
    vecY = zeros(n)

    for k in 1:n
        vecY[k] = vecB[k]
        for i in 1:k-1
            vecY[k] -= matL[k,i]*vecY[i]
        end
    end
    for mink in 1:n
        k = n+1-mink
        vecX[k] = vecY[k]*matL[k,k]
        for i in k+1:n
            vecX[k] -= matL[i,k]*vecX[i]
        end
    end
    return vecX
end

"""
ヤコビアン計算

パラメータ a で与えた二次の画像変換について、ターゲット座標との差のベクトル e の a に関するヤコビアンを計算し、戻り値として返します。
"""
function getYacobian(imgSize = 512, gridSize = 64)
    x = collect(gridSize+0.5:gridSize:imgSize)
    y = collect(gridSize+0.5:gridSize:imgSize)
    n = size(x)[1]
    na = 12

    yacob = Array{Float64}(undef,2*n*n,na)
    for j in 1:n
        for i in 1:n
            idx = i + (j-1)*n
            yacob[2idx-1,1] = -1.0    
            yacob[2idx,1] = 0.0

            yacob[2idx-1,2] = -x[i]    
            yacob[2idx,2] = 0.0
        
            yacob[2idx-1,3] = -y[j]    
            yacob[2idx,3] = 0.0

            yacob[2idx-1,4] = -x[i]^2    
            yacob[2idx,4] = 0.0

            yacob[2idx-1,5] = -x[i]*y[j]    
            yacob[2idx,5] = 0.0

            yacob[2idx-1,6] = -y[j]^2    
            yacob[2idx,6] = 0.0

            yacob[2idx-1,7] = 0.0
            yacob[2idx,7] = -1.0
        
            yacob[2idx-1,8] = 0.0
            yacob[2idx,8] = -x[i]
        
            yacob[2idx-1,9] = 0.0
            yacob[2idx,9] = -y[j]
        
            yacob[2idx-1,10] = 0.0
            yacob[2idx,10] = -x[i]^2
        
            yacob[2idx-1,11] = 0.0
            yacob[2idx,11] = -x[i]*y[j]
        
            yacob[2idx-1,12] = 0.0
            yacob[2idx,12] = -y[j]^2
        end
    end
    return yacob
end

"""
エラーベクトルの計算

1枚目の画像で設定する `gridx` `gridy` から2枚目の画像の `targetX` `targetY` を目的に変換して得た `procX` `procY` を計算し、`targetX` `targetY` と `procX` `procY` の差を計算したベクトル `errorVec` を返します。ヤコビアン取得のため `gridx` `gridy` も同時に返します。
"""
function getErrorVec(vecMap, coefa, gridSize = 64, imgSize = 512)
    n = div(imgSize,gridSize)-1
    gridx = collect(gridSize+0.5:gridSize:imgSize)
    gridy = collect(gridSize+0.5:gridSize:imgSize)
    targetX = Array{Float64}(undef,n*n)
    targetY = Array{Float64}(undef,n*n)
    procX = Array{Float64}(undef,n*n)
    procY = Array{Float64}(undef,n*n)

    for y in 1:n
        for x in 1:n
            targetX[x + n*(y-1)] = gridx[x] + vecMap[y,x,1]
            targetY[x + n*(y-1)] = gridy[y] + vecMap[y,x,2]
            procX[x + n*(y-1)] = coefa[1] + coefa[2]*gridx[x] + coefa[3]*gridy[y] + coefa[4]*gridx[x]^2 + coefa[5]*gridx[x]*gridy[y] + coefa[6]*gridy[y]^2
            procY[x + n*(y-1)] = coefa[7] + coefa[8]*gridx[x] + coefa[9]*gridy[y] + coefa[10]*gridx[x]^2 + coefa[11]*gridx[x]*gridy[y] + coefa[12]*gridy[y]^2
        end
    end

    errorVec = Array{Float64}(undef,2*n*n)

    for idx in 1:n*n
        errorVec[2*idx-1] = targetX[idx] - procX[idx]
        errorVec[2*idx] = targetY[idx] - procY[idx]
    end

    return errorVec
end

"""
画像配列 `img` を、変換係数 `coefa` により画像変換し、返します。
"""
function getNewImage(img, coefa)
    n = size(img)[1]
    bkg = mean(img)
    refX = Array{Int}(undef,n*n)
    refY = Array{Int}(undef,n*n)
    out = Array{Float64}(undef,n,n)
    
    for i in 1:n
        for j in 1:n
            refX[(i-1)*n+j] = Int(round(coefa[1] + coefa[2]*j + coefa[3]*i + coefa[4]*j^2 + coefa[5]*i*j + coefa[6]*i^2))
            refY[(i-1)*n+j] = Int(round(coefa[7] + coefa[8]*j + coefa[9]*i + coefa[10]*j^2 + coefa[11]*i*j + coefa[12]*i^2))
        end
    end

    for i in 1:n
        for j in 1:n
            if (refX[(i-1)*n+j]>=1) && (refX[(i-1)*n+j]<=n) && (refY[(i-1)*n+j]>=1) && (refY[(i-1)*n+j]<=n)
                out[i,j] = img[refY[(i-1)*n+j],refX[(i-1)*n+j]]
            else
                out[i,j] = bkg
            end
        end
    end
    return out 
end

"""
`path` の画像をモノクロ画像として読み込み、画像配列を返します。

画像は 0~1 の Float64 行列として扱われます。
"""
function loadholo(path)
    out = Float32.(channelview(Gray.(load(path))))
end



function CuGetVector!(vecArray::CuDeviceArray{Float32,3},corArray::CuDeviceArray{Float32,2},gridNum::Int64,corArrSize::Int64,intrSize::Int64)
    gridIdxx = (blockIdx().x-1)*blockDim().x + threadIdx().x
    gridIdxy = (blockIdx().y-1)*blockDim().y + threadIdx().y


    if gridIdxx <= gridNum-1 && gridIdxy <= gridNum-1
        x0::Int64 = 0
        y0::Int64 = 0
        
        tmp::Float32 = 0.0
        for i in 1:corArrSize
            for j in 1:corArrSize
                if corArray[corArrSize*(gridIdxy-1)+i,corArrSize*(gridIdxx-1)+j] > tmp
                    x0 = corArrSize*(gridIdxx-1)+j
                    y0 = corArrSize*(gridIdxy-1)+i
                    tmp = corArray[corArrSize*(gridIdxy-1)+i,corArrSize*(gridIdxx-1)+j]
                end
            end
        end

        valy1x0::Float32 = corArray[y0+1,x0]
        valy0x0::Float32 = corArray[y0,x0]
        valyInv1x0::Float32 = corArray[y0-1,x0]
        valy0x1::Float32 = corArray[y0,x0+1]
        valy0xInv1::Float32 = corArray[y0,x0-1]

        if (valy1x0-2.0*valy0x0+valyInv1x0 == 0.0) || (valy0x1-2.0*valy0x0+valy0xInv1 == 0.0)
            valy0x0 += 0.00001
        end

        vecArray[gridIdxy, gridIdxx,1] = Float32(x0) - (valy0x1 - valy0xInv1)/(valy0x1-2.0*valy0x0+valy0xInv1)/2.0 - Float32(intrSize)/2.0 -1.0  - (gridIdxx-1)*corArrSize
        vecArray[gridIdxy, gridIdxx,2] = Float32(y0) - (valy1x0 - valyInv1x0)/(valy1x0-2.0*valy0x0+valyInv1x0)/2.0 - intrSize/2.0 -1.0  - (gridIdxy-1)*corArrSize
    end
    return nothing
end

function CuGetCrossCor!(corArray::CuDeviceArray{Float32,2},img1::CuDeviceArray{Float32,2},img2::CuDeviceArray{Float32,2},gridIdxy::Int64,gridNum::Int64,srchSize::Int64,intrSize::Int64,gridSize::Int64)
    x = (blockIdx().x-1)*blockDim().x + threadIdx().x
    y = (blockIdx().y-1)*blockDim().y + threadIdx().y


    if x <= (srchSize-intrSize+1)*(gridNum-1) && y <= (srchSize-intrSize+1)
        gridIdxx = div(x-1,(srchSize-intrSize+1))+1
        idxx = x - (gridIdxx-1)*(srchSize-intrSize+1)
        idxy = y

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

        corArray[y+(srchSize-intrSize+1)*(gridIdxy-1),x] = num/(CUDA.sqrt(denomA)*CUDA.sqrt(denomB))
    end
    return nothing
end

"""
GPU CUDA Function
画像配列 `image1` と `image2` のPIVマップを計算し、返します。1はx, 2はy。
"""
function getPIVMap_GPU(image1, image2, imgLen = 512, gridSize = 64, intrSize = 64, srchSize = 128)
    gridNum = div(imgLen, gridSize)
    corArray = CuArray{Float32}(undef,((srchSize-intrSize+1)*(gridNum-1),(srchSize-intrSize+1)*(gridNum-1)))
    vecArray = CuArray{Float32}(undef,(gridNum-1,gridNum-1,2))

    blockSize = 16
    threads1 = (blockSize,blockSize)
    blocks1 = (cld((srchSize-intrSize+1)*(gridNum-1),blockSize),cld((srchSize-intrSize+1),blockSize))
    blocks2 = (cld(gridNum-1,blockSize),cld(gridNum-1,blockSize))

    d_img1 = cu(image1)
    d_img2 = cu(image2)

    for idx in 1:gridNum-1
        @cuda threads = threads1 blocks = blocks1 CuGetCrossCor!(corArray,d_img1,d_img2,idx,gridNum,srchSize,intrSize,gridSize)
        synchronize()
        # sleep(0.1)
    end

    # display(corArray)

    @cuda threads = threads1 blocks = blocks2 CuGetVector!(vecArray,corArray,gridNum,srchSize-intrSize+1,intrSize)
    synchronize()
    output = Array(vecArray)
    
    return output
end

function main()
    img1 = loadholo("../PIV/imposed1.jpg")
    # img1 = loadholo("./cam1.bmp")
    # println(img1[1,1])
    img2 = loadholo("../PIV/imposed2.jpg")
    # img2 = loadholo("./cam2.bmp")
    # img1 = img1[1:512,1:512]
    # img2 = img2[1:512,1:512]
    # img1 = loadholo("./cropimg1.png")
    # img2 = loadholo("./cropimg2.png")
    # img1 = rotr90(loadholo("./cropimg1.png"))
    # img2 = rotr90(loadholo("./cropimg2.png"))
    imgLen = size(img1)[1]
    n = 7
    vecArray = getPIVMap_GPU(img1,img2,imgLen);
    # vecArray = getPIVMap_GPU(img1,img2,imgLen);
    # display(vecArray)

    f = Figure(resolution = (1600,500),figure_padding = 1)
    arrowax = Makie.Axis(f[1,3], aspect = 1 , yreversed=false, backgroundcolor="white")
    imageax1 = Makie.Axis(f[1, 1], aspect = DataAspect(), yreversed = false, title = "Camera 1")
    imageax2 = Makie.Axis(f[1, 2], aspect = DataAspect(), yreversed = false, title = "Camera 2")
    imgObservable1 = Observable(rotr90(RGB.(img1,img1,img1)))
    imgObservable2 = Observable(rotr90(RGB.(img2,img2,img2)))
    image!(imageax1,imgObservable1)
    image!(imageax2,imgObservable2)

    vecxObservable = Observable(rotr90(vecArray[:,:,1]))
    vecyObservable = Observable(-rotr90(vecArray[:,:,2]))
    strObservable = Observable(vec(sqrt.(vecArray[:,:,1].^2 .+ vecArray[:,:,2].^2)))

    xs = [i*64 for i in 1:n]
    ys = [i*64 for i in 1:n]
    arrows!(arrowax, xs,ys, vecxObservable,vecyObservable, arrowsize=10, lengthscale=20, arrowcolor = strObservable, linecolor = strObservable)
    # display(f)

    Makie.save("./figure1.pdf",f)

    coefa = fill(1.0,12)
    itr = 1
    yacobian = getYacobian(imgLen)
    hMat = transpose(yacobian)*yacobian

    while itr <= 10
        println("Iteration: ",itr)
        errorVec= getErrorVec(vecArray,coefa)
        display(errorVec)
        println(mean(errorVec.*errorVec))
        # vecxObservable[],vecyObservable[],strObservable[] = errorVecReshape(errorVec,n)
        println()
        # display(hMat)
        deltaCoefa = simEqSolver(choleskyDecom(hMat),yacobian,errorVec)
        # display(deltaCoefa)
        coefa += deltaCoefa
        itr += 1
        sleep(0.1)
    end
    # display(coefa)
    
    img3 = getNewImage(img2, coefa)
    sleep(0.1)
    vecArray = getPIVMap_GPU(img1,img3,imgLen)
    vecxObservable[] = rotr90(vecArray[:,:,1])
    vecyObservable[] = -rotr90(vecArray[:,:,2])
    strObservable[] = vec(sqrt.(vecArray[:,:,1].^2 .+ vecArray[:,:,2].^2))
    imgObservable2[] = rotr90(RGB.(img3,img3,img3))
    Makie.save("./figure1.pdf",f)

    Images.save("./output.png",img3)
end

# main()