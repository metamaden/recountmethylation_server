#!/usr/bin/env python3

""" start_server.py

    Authors: Sean Maden, Abhi Nellore

    Code to start processes used by the server.

"""

import sys
import os
sys.path.insert(0, os.path.join("recountmethylation_server","src"))
import settings
settings.init()

def start_mongodb():
    """
        Start the MongoDB process
    """
    try:
        subprocess.call(settings.runmongocmd.split(" "), shell=False)
        return True
    except:
        return False

def start_rabbitmq():
    """
        Start the RabbitMQ process
    """
    try:
        subprocess.call(settings.runrabbitmqcmd.split(" "), shell=False)
        return True
    except:
        return False

def start_celery():
    """
        Start the celery process
    """
    try:
        subprocess.call(settings.runcelerycmd.split(" "), shell=False)
        return True
    except:
        return False

if __name__ == "__main__":
    """
        Start processes used by the server.
    """
    statuslist = []
    try:
        statuslist.append(start_mongodb())
        statuslist.append(start_rabbitmq())
        statuslist.append(start_start_celery())
        if statuslist[0]==statuslist[0]==statuslist[0]==True:
            try:
                input("All Recount Methylation Server resources were "
                    +"successfully launched! Press ENTER to launch server.py "
                    +"now.")
                subprocess.call(settings.launchserverpycmd.split(" "), 
                        shell=False
                    )
            except SyntaxError:
                print("Finished launching resources for Recount Methylation "
                        +"Server."
                    )
                pass
        else:
            print("Failed to launch recount methylation resources.")
    except:
        print("Failed to launch recount methylation resources.")
