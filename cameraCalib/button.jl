using Gtk, Gtk.ShortNames

# 画面レイアウトの作成
win = Window("Count Click")
v = Box(:v)
l = Label("")
b = Button("Start")
l2 = Label("")
b2 = Button("Stop")
l3 = Label("")
push!(win, v)
push!(v,l)
push!(v, b)
push!(v,l2)
push!(v, b2)
push!(v,l3)
# set_gtk_property!(v, :expand, b, true)

showall(win)

function start()
    println("start")
end

function stop()
    println("stop")
end

# ボタンを押したときに呼ばれる関数
nclick = 0
function click()
    global nclick += 1
    test()
    return nothing
end

# ボタンと関数をつなぐ
signal_connect(x -> start(), b, "clicked")
signal_connect(x -> stop(), b2, "clicked")


if !isinteractive()
    c = Condition()
    signal_connect(win, :destroy) do widget
        notify(c)
    end
    wait(c)
end