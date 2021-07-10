import telebot
from telebot import types
import parsing
from comp_titles import comp_titles

# token: 1728747956:AAEaPlGiZ8gtS_eJrfp0gcXM4LU3NqINhQc

TOKEN = '1707158930:AAEtd_UVsc7La5EeUBXc69VH8fmlwk2e92U'

bot = telebot.TeleBot(TOKEN)
anime = {}
season = 0
episode = 0
seasons_count = 0
episodes_count = 0


def get_link() -> str:
	link = parsing.URL
	if anime is not None:
		link += anime['link'][1:]
	if season != 0:
		link += parsing.SEASON + str(season)
	if episode != 0:
		link += parsing.EPISODE + str(episode) + '.html'
	return link


@bot.message_handler(commands=['start'])
def send_welcome(message):
	# bot.send_message(message.chat.id, 'Подгружаю всю нужную информацию\nсекунду...')
	bot.send_message(message.chat.id, 'Привет!\nЗдесь ты можешь получить серии своего любимого аниме без рекламы(:')
	bot.send_message(message.chat.id, 'Напиши мне название аниме, которое хочешь посмотреть')


@bot.message_handler(func=lambda m: True)
def find_title(message):
	global season
	global episode
	season = 0
	episode = 0
	has_found = -1
	global anime
	for title in parsing.titles:
		result = comp_titles(message.text, title['name'])
		if result == 'maybe':
			has_found = 0
			anime = title
		elif result == title['name']:
			has_found = 1
			anime = title
			break

	if has_found == 1:
		title_success(message)
	elif has_found == 0:
		title_maybe_found(message)
	else:
		title_failure(message)


def title_maybe_found(message):
	keyboard = types.InlineKeyboardMarkup()
	key_yes = types.InlineKeyboardButton(text='Да!', callback_data='yes')
	key_no = types.InlineKeyboardButton(text='Нет):', callback_data='no')
	keyboard.add(key_yes)
	keyboard.add(key_no)
	bot.reply_to(message, f"Ты имел в виду '{anime['name']}'?", reply_markup=keyboard)


@bot.callback_query_handler(func=lambda call: True)
def callback_maybe_found(call):
	if call.data == 'yes':
		title_success(call.message)
	elif call.data == 'no':
		title_failure(call.message)


def title_success(message):
	bot.reply_to(message, f"Нашел!\n{get_link()}")
	global seasons_count
	seasons_count = parsing.how_many_seasons(get_link())
	if seasons_count > 0:
		bot.send_message(message.chat.id, 'Напиши нормер сезона')
		bot.register_next_step_handler(message, find_season)
	else:
		bot.send_message(message.chat.id, 'Напиши нормер серии')
		bot.register_next_step_handler(message, find_episode)
		global episodes_count
		episodes_count = parsing.how_many_episodes(get_link())


def title_failure(message):
	bot.reply_to(message, 'Не нашел):\nПопробуй еще раз')


def find_season(message):
	global season
	season = 0
	success = True
	if season == 0:
		try:
			season = int(message.text)
		except ValueError:
			success = False
			bot.send_message(message.chat.id, 'Цифрами(:')
			bot.register_next_step_handler(message, find_season)

	flag = season <= seasons_count and bool(season)
	if flag:
		season_success(message)
	elif success:
		bot.reply_to(message, 'Похоже, этот сезон еще не вышел):\nМожешь посмотреть другой!')
		bot.register_next_step_handler(message, find_season)


def season_success(message):
	bot.reply_to(message, f"Нашел!\n{get_link()}")
	bot.send_message(message.chat.id, 'Напиши номер серии')
	bot.register_next_step_handler(message, find_episode)
	global episodes_count
	episodes_count = parsing.how_many_episodes(get_link())


def find_episode(message):
	global episode
	episode = 0
	success = True
	if episode == 0:
		try:
			episode = int(message.text)
		except ValueError:
			success = False
			bot.send_message(message.chat.id, 'Цифрами(:')
			bot.register_next_step_handler(message, find_episode)

	flag = episode <= episodes_count and bool(episode)
	if flag:
		episode_success(message)
	elif success:
		bot.reply_to(message, 'Похоже, эта серия еще не вышла):\nМожешь посмотреть другую!')
		bot.register_next_step_handler(message, find_episode)


def episode_success(message):
	bot.reply_to(message, f"Нашел!\n{get_link()}")
	bot.send_message(message.chat.id, 'Приятного просмотра(:')


if __name__ == '__main__':
	parsing.parse_titles()
	bot.polling()
