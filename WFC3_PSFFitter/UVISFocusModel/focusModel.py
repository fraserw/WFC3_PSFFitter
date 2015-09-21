import pickle
from scipy import interpolate as interp

class focusModel:
    def __init__(self):
        """
        Init as focus=focusModel.focusModel()
        """

        with open('/Users/fraserw/git/WFC3_PSFFitter/WFC3_PSFFitter/UVISFocusModel/focusData.pickle') as han:
            
            [mjd,focus]=pickle.load(han)

        self.focus=interp.interp1d(mjd,focus)

    def __call__(self,x):
        """
        Call as focus(mjd)
        -where mjd is a floating point modified julian date between 2003 and 2013
        """

        return self.focus(x)


if __name__=="__main__":
    with open('FocusAll.txt') as han:
        data=han.readlines()

    mjd=[]
    focus=[]
    for i in range(len(data)):
        s=data[i].split()
        mjd.append(float(s[0]))
        focus.append(float(s[len(s)-1]))

    with open('focusData.pickle','w+') as han:
        pickle.dump([mjd,focus],han)
