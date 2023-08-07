<<<<<<< HEAD
debug_muchang = False
debug_michelle = True
=======
debug_muchang = True
debug_michelle = False
>>>>>>> 90f4e5e4036d8333fd80d8f4f92d284bb3a58f80
debug_abhishek = False 

def debug(name:str, error_message:str): 
    if name == "Muchang" and debug_muchang == True: 
        print(error_message)
    if name == "Michelle" and debug_michelle == True: 
        print(error_message) 
    if name == "Abhishek" and debug_abhishek == True: 
        print(error_message)