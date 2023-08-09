debug_muchang = True
debug_michelle = False
debug_abhishek = False 

def debug(name:str, error_message:str): 
    if name == "Muchang" and debug_muchang == True: 
        print(error_message)
    if name == "Michelle" and debug_michelle == True: 
        print(error_message) 
    if name == "Abhishek" and debug_abhishek == True: 
        print(error_message)