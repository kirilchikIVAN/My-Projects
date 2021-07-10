import logic
import tkinter as tk


def _from_rgb(rgb):
    r, g, b = rgb
    r %= 256
    g %= 256
    b %= 256
    return f'#{r:02x}{g:02x}{b:02x}'


class Clicker:
    LEFT_COL = 2
    MID_COL = 325
    RIGHT_COL = 550
    DELTA = 150
    DEFAULT_BG = 'lightgreen'

    class ButtonCommand:
        def __init__(self, level: int, player: logic.Player, window_: tk.Tk, label: tk.StringVar,
                     main_text: tk.Label, button_left_text: list = 0, clicker_left_text: list = 0):
            self.level = level
            self.window = window_
            self.player = player
            self.label = label
            self.main_text = main_text
            self.button_left_text = button_left_text
            self.clicker_left_text = clicker_left_text

        def click_button(self):
            if self.player.buy_clicker(self.level):
                self.label.set('now ' + str(self.player.clickers[self.level]))
                self.main_text['text'] = 'Money: ' + str(self.player.my_money())

        def click_upgrade(self):
            if self.player.clicker_upgrade(self.level):
                self.label.set('level ' + str(self.player.clickers_stages[self.level]))
                self.button_left_text[self.level]['text'] =\
                    'costs: ' + str(self.player.shop.clicker_upgrades[self.level].cost)
                self.main_text['text'] = 'Money: ' + str(self.player.my_money())
                self.clicker_left_text[self.level]['text'] =\
                    'costs: ' + str(self.player.shop.clickers[self.level].cost) + '\n' + \
                    'gains: ' + str(round(self.player.shop.clickers[self.level].gain, 2)) + '/s'

    def __init__(self, window: tk.Tk):
        self.player = logic.Player()
        self.window = window
        self.window['bg'] = self.DEFAULT_BG
        self.window.geometry('800x600')
        self.clicker_colors = [_from_rgb((255, 255 - i, 255 - i)) for i in range(40, 255) if i % 40 == 0]
        self.upgrade_colors = [_from_rgb((255 - i, 255 - i, 255)) for i in range(40, 255) if i % 40 == 0]
        self.main_label_text = tk.Label(self.window)
        self.auto_gain_text = tk.Label(self.window)
        self.clicker_buttons_text = [tk.StringVar() for i in range(self.player.shop.NUMB_OF_CLICKERS + 1)]
        self.clicker_left_text = [tk.Label() for i in range(self.player.shop.NUMB_OF_CLICKERS + 1)]
        self.clicker_buttons = [tk.Button(self.window, padx="10", pady="10",
                                          bg=self.clicker_colors[i % len(self.clicker_colors)],
                                          textvariable=self.clicker_buttons_text[i],
                                          command=self.ButtonCommand(i,
                                                                     self.player,
                                                                     self.window,
                                                                     self.clicker_buttons_text[i],
                                                                     self.main_label_text).click_button)
                                for i in range(self.player.shop.NUMB_OF_CLICKERS)]

        self.main_button_gain = tk.Label()
        self.main_button = tk.Button(self.window, padx='50', pady='50', text='Click',
                                     command=self.click_main_button, bg='green')

        self.upgrade_buttons_text = [tk.StringVar() for i in range(self.player.shop.NUMB_OF_CLICKERS + 1)]
        self.upgrade_left_text = [tk.Label() for i in range(self.player.shop.NUMB_OF_CLICKERS + 1)]
        self.upgrade_buttons = [tk.Button(self.window, padx="10", pady="10",
                                          bg=self.upgrade_colors[i % len(self.upgrade_colors)],
                                          textvariable=self.upgrade_buttons_text[i],
                                          command=self.ButtonCommand(i,
                                                                     self.player,
                                                                     self.window,
                                                                     self.upgrade_buttons_text[i],
                                                                     self.main_label_text,
                                                                     self.upgrade_left_text,
                                                                     self.clicker_left_text).click_upgrade)
                                for i in range(self.player.shop.NUMB_OF_CLICKERS)]
        self.main_upgrade = tk.Button(self.window, padx='40', pady='30', text='level 1',
                                      command=self.click_main_upgrade, bg=_from_rgb((0, 150, 0)))
        self.main_upgrade_cost = tk.Label()
        self.victory = tk.Button(self.window, padx='40', pady='20', text='100 000 to win',
                                 command=self.victory_func, bg='gold')

    def victory_func(self):
        if self.player.money >= 100000:
            self.player.money -= 100000
            self.update_money()
            victory_window = tk.Tk()
            victory_window.geometry('140x30')
            victory_window['bg'] = 'gold'
            victory_label = tk.Label(victory_window, text='!!!You won!!!', bg='gold')
            victory_label.place(x=0, y=0)
            victory_window.mainloop()

    def click_main_upgrade(self):
        if self.player.player_upgrade():
            self.main_upgrade['text'] = 'level ' + str(self.player.my_stage)
            self.main_upgrade_cost['text'] = 'costs: ' + str(self.player.shop.player_upgrade.cost)
            self.main_button_gain['text'] = 'gains: ' + str(int(self.player.self_gain)) + '/click'
            self.update_money()

    def click_main_button(self):
        self.player.make_money()
        self.update_money()

    def init_main_label(self):
        self.main_label_text['bg'] = self.DEFAULT_BG
        self.main_label_text['text'] = 'Money: ' + str(self.player.my_money())
        self.main_label_text.place(x=self.MID_COL, y=50)
        self.auto_gain_text['bg'] = self.DEFAULT_BG
        self.auto_gain_text['text'] = 'Autogain: 0/s'
        self.auto_gain_text.place(x=self.MID_COL, y=80)

    def init_clicker_buttons(self):
        for i in range(self.player.shop.NUMB_OF_CLICKERS):
            self.clicker_buttons_text[i].set('now ' + str(0))
            self.clicker_buttons[i].place(x=self.LEFT_COL+self.DELTA, y=100 + 60 * (i + 1))
            self.clicker_left_text[i]['bg'] = self.DEFAULT_BG
            self.clicker_left_text[i]['text'] = 'costs: ' + str(self.player.shop.clickers[i].cost) + '\n' +\
                                                'gains: ' + str(round(self.player.shop.clickers[i].gain, 2)) + '/s'
            self.clicker_left_text[i].place(x=self.LEFT_COL, y=105 + 60 * (i + 1))

    def init_clicker_upgrades(self):
        for i in range(self.player.shop.NUMB_OF_CLICKERS):
            self.upgrade_buttons_text[i].set('level ' + str(1))
            self.upgrade_buttons[i].place(x=self.RIGHT_COL+self.DELTA, y=100 + 60 * (i + 1))
            self.upgrade_left_text[i]['bg'] = self.DEFAULT_BG
            self.upgrade_left_text[i]['text'] = 'costs: ' + str(self.player.shop.clicker_upgrades[i].cost) + '\n'
            self.upgrade_left_text[i].place(x=self.RIGHT_COL, y=110 + 60 * (i + 1))
        self.main_upgrade.place(x=self.MID_COL, y=410)
        self.main_upgrade_cost['bg'] = self.DEFAULT_BG
        self.main_upgrade_cost['text'] = 'costs: ' + str(self.player.shop.player_upgrade.cost)
        self.main_upgrade_cost.place(x=self.MID_COL, y=380)

    def init_main_button(self):
        self.main_button.place(x=self.MID_COL, y=200)
        self.main_button_gain['bg'] = self.DEFAULT_BG
        self.main_button_gain['text'] = 'gains: ' + str(int(self.player.self_gain)) + '/click'
        self.main_button_gain.place(x=self.MID_COL, y=170)

    def update_money(self):
        self.main_label_text['text'] = 'Money: ' + str(self.player.my_money())
        self.auto_gain_text['text'] = 'Autogain: ' + str(self.player.get_money()) + '/s'

    def update(self):
        self.player.money += self.player.get_money()
        self.update_money()
        self.window.after(1000, self.update)

    def run(self):
        self.init_main_label()
        self.init_clicker_buttons()
        self.init_clicker_upgrades()
        self.init_main_button()
        self.update()
        self.victory.place(x=70, y=500)

        text_auto_clickers = tk.Label(self.window, text='autoclickers:', bg=self.DEFAULT_BG)
        text_auto_clickers.place(x=self.LEFT_COL, y=100)
        text_upgrades = tk.Label(self.window, text="autoclicker's upgrades:", bg=self.DEFAULT_BG)
        text_upgrades.place(x=self.RIGHT_COL, y=100)
        text_click_upgrade = tk.Label(self.window, text="click upgrade:", bg=self.DEFAULT_BG)
        text_click_upgrade.place(x=self.MID_COL, y=350)


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.clicker = Clicker(self)

    def run(self):
        self.clicker.run()
        self.mainloop()


if __name__ == "__main__":
    app = App()
    app.run()
