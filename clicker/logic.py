class Clicker:
    DEFAULT_COST = 50
    DEFAULT_GAIN = 1
    PERCENT = 400

    def __init__(self, level: int = 0):
        self.cost = self.DEFAULT_COST
        self.gain = self.DEFAULT_GAIN
        self.upgrade(level)

    def upgrade(self, times: int):
        self.cost *= 3 ** times
        self.gain *= (self.PERCENT / 100) ** times


class Upgrade:
    DEFAULT_COST = 100
    CLICKER_PERCENTAGE = 30
    PLAYER_PERCENTAGE = 100

    def __init__(self, clicker: bool, stage: int, level = 0):
        self.cost = 0
        self.percentage = 0
        if clicker:
            self.make_clicker_upgrade(level, stage)
        else:
            self.make_player_upgrade(stage)

    def make_clicker_upgrade(self, level: int, stage: int):
        self.cost = self.DEFAULT_COST * 2 ** (level + 1) * stage
        self.percentage = self.CLICKER_PERCENTAGE * stage

    def make_player_upgrade(self, stage: int):
        self.cost = self.DEFAULT_COST * 2 ** (stage - 1)
        self.percentage = self.PLAYER_PERCENTAGE


class Shop:
    NUMB_OF_CLICKERS = 5

    def __init__(self):
        self.clickers = [Clicker(i) for i in range(self.NUMB_OF_CLICKERS)]
        self.clicker_upgrades = [Upgrade(clicker=True, stage=1, level=i) for i in range(self.NUMB_OF_CLICKERS)]
        self.player_upgrade = Upgrade(clicker=False, stage=1)

    def update_clicker_upgrade(self, level: int, new_stage: int):
        self.clicker_upgrades[level].make_clicker_upgrade(level, new_stage)

    def update_player_upgrade(self, new_stage: int):
        self.player_upgrade.make_player_upgrade(new_stage)

    def get_clicker(self, level: int) -> Clicker:
        return self.clickers[level]


class Player:
    DEFAULT_SELF_GAIN = 1
    DEFAULT_GAIN = 0

    def __init__(self):
        self.money = 0
        self.my_stage = 1
        self.self_gain = self.DEFAULT_SELF_GAIN
        self.gain = self.DEFAULT_GAIN
        # self.clicker_max_level = 0
        self.clickers = [0] * Shop.NUMB_OF_CLICKERS
        self.clickers_stages = [1] * Shop.NUMB_OF_CLICKERS
        self.shop = Shop()

    def my_money(self) -> int:
        return int(self.money)

    def enough_money(self, need: int) -> bool:
        return need <= self.money

    def buy_clicker(self, level: int) -> bool:
        new_clicker = self.shop.get_clicker(level)

        if not self.enough_money(new_clicker.cost):
            return False
        self.money -= new_clicker.cost
        if len(self.clickers) == level:
            self.clickers.append(0)
            self.clickers_stages.append(0)

        self.clickers[level] += 1
        return True

    def clicker_upgrade(self, level: int) -> bool:
        new_upgrade = self.shop.clicker_upgrades[level]

        if not self.enough_money(new_upgrade.cost):
            return False

        self.money -= new_upgrade.cost
        self.shop.clickers[level].gain *= 1 + new_upgrade.percentage / 100
        self.clickers_stages[level] += 1
        self.shop.update_clicker_upgrade(level, self.clickers_stages[level])
        return True

    def player_upgrade(self):
        new_upgrade = self.shop.player_upgrade

        if not self.enough_money(new_upgrade.cost):
            return False

        self.money -= new_upgrade.cost
        self.self_gain *= 1 + new_upgrade.percentage / 100
        self.my_stage += 1
        self.shop.update_player_upgrade(self.my_stage)
        return True

    def get_money(self) -> int:
        all_gain = 0
        for i in range(len(self.clickers)):
            all_gain += self.clickers[i] * self.shop.clickers[i].gain
        return int(all_gain)

    def make_money(self):
        self.money += self.self_gain
