using Gtk, Gtk.ShortNames

# 画面レイアウトの作成
# win = Window("Count Click")
# v = Box(:v)
# l = Label("")
# b = Button("Stop")
# l2 = Label("")
# push!(win, v)
# push!(v,l)
# push!(v, b)
# push!(v,l2)

# showall(win)

# function stop()
#     println("stop")
# end

# ボタンを押したときに呼ばれる関数
# state = "running"
# function show(state)
#     while state == "running"
#         println("running",rand())
#         sleep(0.1)
#     end
# end

# function stop(state)
#     global state = "stop"
# end

# ボタンと関数をつなぐ
# signal_connect(x -> stop(state), b, "clicked")

# show(state)

function main()
    win = Window("Count Click")
    v = Box(:v)
    l = Label("")
    b = Button("Stop")
    l2 = Label("")
    push!(win, v)
    push!(v,l)
    push!(v, b)
    push!(v,l2)

    showall(win)


    state = "running"
    while state == "running"
        println("running",rand())
        sleep(0.1)
        # signal_connect(x -> state = "stop", b, "clicked")
    end

    signal_connect(x -> state = "stop", b, "clicked")


    if !isinteractive()
        c = Condition()
        signal_connect(win, :destroy) do widget
            notify(c)
        end
        wait(c)
    end
end

main()