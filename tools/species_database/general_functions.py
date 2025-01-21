import re

class myColors:
    """ Add different colors and styles (ANSI code) to strings """
    endc = '\x1b[0m'
    black = '\x1b[30m'
    red = '\x1b[31m'
    green = '\x1b[32m'
    yellow = '\x1b[33m'
    blue = '\x1b[34m'
    purple = '\x1b[35m'
    cyan = '\x1b[36m'
    white = '\x1b[37m'
    bold = '\x1b[1m'
    underlinE = '\x1b[4m'

def bold(text):
    return myColors.bold + text + myColors.endc

def underlinE(text):
    return myColors.underlinE + text + myColors.endc

def red(text):
    return myColors.red + text + myColors.endc

def green(text):
    return myColors.green + text + myColors.endc

def blue(text):
    return myColors.blue + text + myColors.endc

def yellow(text):
    return myColors.yellow + text + myColors.endc

def purple(text):
    return myColors.purple + text + myColors.endc

def cyan(text):
    return myColors.cyan + text + myColors.endc

def get_valid_input(prompt, validator_func):
    try:
        while True:
            user_input = input(prompt)
            if validator_func(user_input):
                return user_input
            else:
                print(bold(red('Invalid input!')))
    except KeyboardInterrupt:
        print('\nKeyboardInterrupt: Program terminated by user.')

def create_prompt(*args):
    prompt_for_user = bold('\nPlease enter') + '\n '
    for (i, arg) in enumerate(args):
        prompt_for_user = prompt_for_user + purple(str(i + 1)) + ' ' + arg + ' ' + ' or \n '
    prompt_for_user = prompt_for_user + purple(str(i + 2)) + ' to exit program\n-->'
    return prompt_for_user

def own_exit():
    print(bold(red('Exiting')))
    exit(1)

def int_to_Roman(num):
    val = (1000, 900, 500, 400, 100, 90, 50, 40, 10, 9, 5, 4, 1)
    syb = ('M', 'CM', 'D', 'CD', 'C', 'XC', 'L', 'XL', 'X', 'IX', 'V', 'IV', 'I')
    roman_num = ''
    for i in range(len(val)):
        count = int(num / val[i])
        roman_num += syb[i] * count
        num -= val[i] * count
    return roman_num

def remove_from_list(List, *args_to_remove):
    NewList = []
    for element in List:
        if element == 'electron':
            NewList.append('el')
        elif element not in args_to_remove:
            NewList.append(element)
    return NewList

def get_interaction_id(species_name):
    if sum((1 for c in species_name.replace('Ion', '') if c.isupper())) == 1 and (not bool(re.search('\\d+', re.sub('Ion\\d+', '', species_name)))):
        if not bool(re.search('[A-Za-z]*\\d', re.sub('Ion\\d+', '', species_name))):
            if 'Ion' in species_name:
                interactionID = 10
            elif 'Ion' not in species_name:
                interactionID = 1
    elif bool(re.search('\\d+', re.sub('Ion\\d+', '', species_name))) or sum((1 for c in species_name.replace('Ion', '') if c.isupper())) != 1:
        if 'Ion' in species_name:
            interactionID = 20
        elif 'Ion' not in species_name:
            interactionID = 2
    elif 'electron' in species_name:
        interactionID = 4
    return interactionID