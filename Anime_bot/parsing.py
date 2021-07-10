import requests
from bs4 import BeautifulSoup
import re


URL = 'https://jut.su/'
TITLE_URL = 'https://jut.su/anime/'
SEASON = 'season-'
EPISODE = '/episode-'
HEADERS = {
    'user-agent': 'Mozilla/5.0 (X11; Ubuntu; Linux x86_64; rv:88.0) Gecko/20100101 Firefox/88.0',
    'accept': '*/*',
}
MAX_PAGE = 26
# 27

titles = []


def get_html(url, params=None):
    r = requests.get(url, headers=HEADERS, params=params)
    return r


def get_content(url, params=None):
    html = requests.get(url, headers=HEADERS, params=params).text
    return BeautifulSoup(html, 'html.parser')


def get_titles(html):
    soup = BeautifulSoup(html, 'html.parser')
    items = soup.find_all('div', class_='all_anime')
    for i in range(len(items)):
        items[i] = items[i].previous_element

    global titles
    for item in items:
        titles.append({
            'name': item.find('div', class_='aaname').text,
            'link': item.get('href')
        })


def parse_titles():
    for i in range(1, MAX_PAGE + 1):
        print(f'Parsing page №{i}...')
        html = get_html(TITLE_URL + f'page-{i}')
        if html.status_code == 200:
            get_titles(html.text)
        else:
            print('Error')


def how_many_seasons(link) -> int:
    soup = get_content(link)
    count = 0
    for item in soup.find_all('h2', class_="the-anime-season"):
        if len(re.findall(r'.*сезон', item.text)) == 1:
            count += 1
    return count


def how_many_episodes(link) -> int:
    soup = get_content(link)
    return len(soup.find_all('a', class_='the_hildi'))


if __name__ == '__main__':
    # parse_series('/shingeki-kyojin/', '10')
    pass
