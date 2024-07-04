###################################################################################################
#   General functions
###################################################################################################
class myColors :
    """ Add different colors and styles (ANSI code) to strings """
    # After coloring, change back to \033
    endc   = '\033[0m'

    # Regular Colors
    black  = '\033[30m'
    red    = '\033[31m'
    green  = '\033[32m'
    yellow = '\033[33m'
    blue   = '\033[34m'
    purple = '\033[35m'
    cyan   = '\033[36m'
    white  = '\033[37m'

    # Text Style
    bold = '\033[1m'
    underlinE = '\033[4m'

def bold(text) :
    return myColors.bold+text+myColors.endc
def underlinE(text) :
    return myColors.underlinE+text+myColors.endc
def red(text) :
    return myColors.red+text+myColors.endc
def green(text) :
    return myColors.green+text+myColors.endc
def blue(text) :
    return myColors.blue+text+myColors.endc
def yellow(text) :
    return myColors.yellow+text+myColors.endc
def purple(text) :
    return myColors.purple+text+myColors.endc
def cyan(text) :
    return myColors.cyan+text+myColors.endc

# input function
def get_valid_input(prompt, validator_func):
    try:
        while True:
            user_input = input(prompt)
            if validator_func(user_input):
                return user_input
            else:
                print(bold(red("Invalid input!")))
    except KeyboardInterrupt:
        print("\nKeyboardInterrupt: Program terminated by user.")

# create prompt for input functions
def create_prompt(*args):
    prompt_for_user = bold('\nPlease enter') + "\n "
    for i,arg in enumerate(args):
        prompt_for_user = prompt_for_user + purple(str(i+1)) + " " + arg + " " +" or \n "
    prompt_for_user = prompt_for_user + purple(str(i+2)) + " to exit program\n-->"
    return prompt_for_user

def own_exit():
    print(bold(red("Exiting")))
    exit(1)

# Function to determine the Roman number from an integer
def int_to_Roman(num):
   val = (1000, 900,  500, 400, 100,  90, 50,  40, 10,  9,   5,  4,   1)
   syb = ('M',  'CM', 'D', 'CD','C', 'XC','L','XL','X','IX','V','IV','I')
   roman_num = ""
   for i in range(len(val)):
      count = int(num / val[i])
      roman_num += syb[i] * count
      num -= val[i] * count
   return roman_num

def remove_from_list(List, *args_to_remove):
    NewList = []
    for element in List:
        if element == 'electron':   # special case for reaction handling
            NewList.append('el')
        elif element not in args_to_remove:
            NewList.append(element)
    return NewList

def parse_selection(selection, total_files):
    selected_indices = set()
    parts = selection.split(',')
    for part in parts:
        if '-' in part:
            start, end = part.split('-')
            start, end = int(start.strip()), int(end.strip())
            selected_indices.update(range(start, end + 1))
        else:
            selected_indices.add(int(part.strip()))
    return [idx - 1 for idx in selected_indices if 0 <= idx - 1 < total_files]

def read_file_by_numbers(txt_files):
    while True:
        try:
            file_numbers = input("Enter the numbers corresponding to the files you want to read (comma separated, ranges allowed): ")
            selected_indices = parse_selection(file_numbers, len(txt_files))
            if selected_indices:
                selected_files = [txt_files[i] for i in selected_indices]
                break
            else:
                print("Invalid numbers. Please try again.")
        except ValueError:
            print("Invalid input. Please enter numbers separated by commas or ranges.")
    return selected_files