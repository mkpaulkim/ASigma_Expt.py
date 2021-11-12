import hardware as hw

uno = hw.Arduino()
LED = 3     # pin 11

# switch = 'ON'

while True:
    switch = input('on or off: ')

    if switch.lower() == 'on':
        uno.set_LED(LED)
    else:
        uno.set_LED(0)







