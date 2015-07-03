import numpy as np
import cv2
"""
cap = cv2.VideoCapture('/home/andrew/pupil/recordings/2015_06_23/001/eye0.mkv')

a = True
while(a == True):
    # Capture frame-by-frame
    ret, frame = cap.read()

    # Our operations on the frame come here
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

    # Display the resulting frame
    cv2.imshow('frame',gray)
    if cv2.waitKey(1) & 0xFF == ord('q'):
        a = False

# When everything done, release the capture
cap.release()
cv2.destroyAllWindows()
"""


##################


cap = cv2.VideoCapture("/home/andrew/pupil/recordings/2015_06_23/001/eye0.mkv")
idx = 0
while True:
    flag, frame = cap.read()
    if flag:
        # The frame is ready and already captured
        cv2.imshow('video', frame)
        cv2.waitKey(1)
        cv2.imwrite('andrew_eye_%000i.png'%idx,frame)
        idx +=1
    else:
    	break
    