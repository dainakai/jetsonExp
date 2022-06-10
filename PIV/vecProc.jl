using Glob
using DelimitedFiles

function main()
    date = ARGS[1]
    vecPath = "./vectorResult/"*date*"/*/*.dat"
    vecArrPath = glob(vecPath)
    vecArray = readdlm(vecArrPath[1])
    for i in 2:length(vecArrPath)
        tmp = readdlm(vecArrPath[i])
        for idx in 1:7^2
            if tmp[idx,3]^2+tmp[idx,4]^2 < vecArray[idx,3]^2+vecArray[idx,4]^2
                vecArray[idx,3] = tmp[idx,3]
                vecArray[idx,4] = tmp[idx,4]
            end
        end
    end
    display(vecArray)
    println("")
    writedlm("./vecArray.dat",vecArray)
end

main()