def comp_titles(message_: str, title_: str):
    message = message_.lower()
    title = title_.lower()
    different_len = len(message) != len(title)

    dif1, dif2 = 0, 0
    len_ = min(len(message), len(title))
    for i in range(len_):
        if message[i] != title[i]:
            dif1 += 1
        if message[len(message) - 1 - i] != title[len(title) - 1 - i]:
            dif2 += 1

    if dif1 / len_ > 0.5 and dif2 / len_ > 0.5:
        return 'not correct'
    elif different_len:
        return 'maybe'
    return title_
