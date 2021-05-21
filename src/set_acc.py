#!/usr/bin/env python3

""" set_acc.py

    Author: Sean Maden

    Description:
    Set the platform accession ID for the recountmethylation instance. Contains 
    user prompts to specify the target accession ID for the instance.

"""

import os, sys
sys.path.insert(0, os.path.join("recountmethylation_server","src"))

def acc_dialogue(current_platform):
    """ acc_dialogue
    
        Arguments:
        * current_platform: Current platform accession ID specified in 
            settings.py

        Returns:
        * User-specified platform accession ID, based on prompt outcomes.

    """
    input2 = current_platform
    prompt1 = "".join(["(Y/N) Current target platform accession ID is ", 
        current_platform, ".\n", "Do you want to change the target platform?"])
    prompt2 = " ".join(["(1/2/3) Enter the number of the platform to target:\n",
        "1. HM450K (GPL13534)\n", "2. EPIC/HM850K (GPL21145)\n", 
        "3. HM27K (GPL8490)\n"])
    input1 = input(prompt1)
    if input1 in ["y", "Y", "yes", "Yes", "YES"]:
        input2 = input(prompt2)
        if not input2 in ["1", "2", "3"]:
            print("Error, invalid input. Canceling prompt.")
        else:
            if input2 == "1":
                input2 = "GPL13534"
            elif input2 == "2":
                input2 = "GPL21145"
            else:
                input2 = "GPL8490"
    elif input1 not in ["n", "N", "no", "No", "NO"]:
        print("Invalid input. Canceling prompt.")
    else:
        print("Canceling prompt...")
    return input2

def get_acc_opt(settings_path):
    """ get_acc_opt

        Get the accession ID options.

        Arguments
        * settings_path: Path to settings.py script.

        Returns 
        * True, sets the accession in settings.py based on prompts.
    """
    with open(settings_path, "r") as rf:
        lines = rf.readlines()
    acc_line = lines[24]
    acc_split = acc_line.split("'")
    acc_split[1] = acc_dialogue(acc_split[1])
    lines[24] = "'".join(acc_split)
    with open(settings_path, "w") as wf:
        wf.writelines(lines)
    return True

if __name__ == "__main__":
    src_path = os.path.join("recountmethylation_server", "src")
    settings_path = os.path.join(src_path,"settings.py")
    get_acc_opt(settings_path)


