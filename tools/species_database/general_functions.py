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