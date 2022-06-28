import time


def String2FloatList(string):
    L = string.split()
    n = len(L)
    for i in range(n):
        L[i] = float(L[i])
    return L


def progressBar(step, limit, deltaT, width=30):
    percent = step / limit * 100
    left = int(width * percent / 100)
    right = width - left
    print('\r[', '#' * left, ' ' * right, ']',
          f' {percent:.0f}%', f'   {step:.0f}|{limit:.0f}',
          f'   {deltaT:.2f} s',
          sep='', end='', flush=True)


def sessionName(separator=''):
    t0 = time.time()
    struct = time.localtime(t0)
    string = str(struct.tm_year) + separator
    nMonths = str(struct.tm_mon)  # MONTHS
    if len(nMonths) == 1:
        nMonths = '0' + nMonths
    string = string + nMonths + separator
    nDays = str(struct.tm_mday)  # DAYS
    if len(nMonths) == 1:
        nDays = '0' + nDays
    string = string + nDays + separator
    nHours = str(struct.tm_hour)  # HOURS
    if len(nHours) == 1:
        nHours = '0' + nHours
    string = string + nHours + separator
    nMinutes = str(struct.tm_min)  # MINUTES
    if len(nMinutes) == 1:
        nMinutes = '0' + nMinutes
    string = string + nMinutes + separator
    nSecs = str(struct.tm_sec)  # SECONDS
    if len(nSecs) == 1:
        nSecs = '0' + nSecs
    string = string + nSecs
    return string
