import wx 
import sys

class SliderFrame(wx.Frame):
    tileOfWindow='Scale'
    titleOfScale='Scale'

    tileOfLeft='Left'
    titleOfRight='Right'

    windowWidth=300

    inputText=[]
    inputText2=[]
    slider=[]
	
    scaleMin=0
    scaleMax=10
	
    numOfIteration=0
    thisIteration=0

    margin=100

    button=[]
    panel=[]

    values=[]
	
    def readLineFromFile(self,f): # removes all unwanted "Enter" charachters and returns real lines
        myLine=chr(10)
        while ord(myLine[0])==10:
            myLine=f.readline()
        return myLine		
	
    def buttonPressed(self,event):
        self.values.append(self.slider.GetValue())
        self.thisIteration=self.thisIteration+1
        if self.thisIteration==self.numOfIteration:
            self.button.SetLabel('Save')	
        if self.thisIteration>self.numOfIteration:
            self.saveAndClose()

        self.numOfIteration=self.readFromFile(self.thisIteration) # on init get data from config.txt

        self.scaleLabel.Destroy()
        self.leftLabel.Destroy()
        self.rightLabel.Destroy()
        self.slider.Destroy()

        panel = wx.Panel(self, -1)

        self.scaleLabel = wx.StaticText(self.panel, -1, self.titleOfScale, (10, 100), (self.windowWidth-20, -1), wx.ALIGN_CENTER)
        self.leftLabel = wx.StaticText(self.panel, -1, self.tileOfLeft,(0,150),(self.margin*2,-1),wx.ALIGN_CENTER)
        self.rightLabel = wx.StaticText(self.panel, -1, self.titleOfRight,(self.windowWidth-self.margin*2,150),(self.margin*2,-1),wx.ALIGN_CENTER)

        self.slider = wx.Slider(self.panel, (self.scaleMax-self.scaleMin)/2, (self.scaleMax-self.scaleMin)/2, self.scaleMin, self.scaleMax, pos=(self.margin, 170),
                size=(self.windowWidth-self.margin*2, -1),
                style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS   )
        self.slider.SetTickFreq((self.scaleMax-self.scaleMin)/2, (self.scaleMax-self.scaleMin)/2)
        	
    def saveAndClose(self): 
        filename=self.inputText.GetValue()+'_'+self.inputText2.GetValue()+'_'+self.inputText3.GetValue()+'.txt' # filename made of subject ID +subject session +measurment
        fw=open(filename,'w')
        for itr in range(0,len(self.values)):
            fw.write(str(self.values[itr])+'\n')
            print (str(self.values[itr])+'\n')
        fw.close()
        sys.exit(0)
		
    def readFromFile(self,iteration): # get initial data from config.txt
            f = open('config.txt', 'r')
            numOfIterations=int(self.readLineFromFile(f))
            for itr in range(1,iteration+1):
	            self.tileOfWindow=self.readLineFromFile(f)  # this line reads the title for the window
	            self.windowWidth=int(self.readLineFromFile(f))  # this line reads the width of the window
	            self.scaleMin=int(self.readLineFromFile(f))  # this line reads the min of the scale
	            self.scaleMax=int(self.readLineFromFile(f))  # this line reads the max of the scale
	            self.titleOfScale=self.readLineFromFile(f)  # this line reads the title for the scale
	            self.tileOfLeft=self.readLineFromFile(f)  # this line reads the title for the left side
	            self.titleOfRight=self.readLineFromFile(f)  # this line reads the title for the right side

            return numOfIterations
			
    def __init__(self):
        self.numOfIteration=self.readFromFile(1) # on init get data from config.txt
        self.thisIteration=1
        wx.Frame.__init__(self, None, -1,self.tileOfWindow, size=(self.windowWidth, 350)) # creates visual window
        self.panel = wx.Panel(self, -1)

		# creates labels + input boxes
        inputLabel = wx.StaticText(self.panel, -1, "Subject ID:",(10,20)) 
        self.inputText = wx.TextCtrl(self.panel, -1, "1", (100,20),size=(175, -1))
        self.inputText.SetInsertionPoint(0)

        inputLabel2 = wx.StaticText(self.panel, -1, "Session:",(10,40))
        self.inputText2 = wx.TextCtrl(self.panel, -1, "1", (100,40),size=(175, -1))
        self.inputText2.SetInsertionPoint(0)

        inputLabel3 = wx.StaticText(self.panel, -1, "Measurement:",(10,60))
        self.inputText3 = wx.TextCtrl(self.panel, -1, "1", (100,60),size=(175, -1))
        self.inputText3.SetInsertionPoint(0)

		# creates scale and labels to match
		
        self.scaleLabel = wx.StaticText(self.panel, -1, self.titleOfScale, (10, 100), (self.windowWidth-20, -1), wx.ALIGN_CENTER)
		
        self.leftLabel = wx.StaticText(self.panel, -1, self.tileOfLeft,(0,150),(self.margin*2,-1),wx.ALIGN_CENTER)
        self.rightLabel = wx.StaticText(self.panel, -1, self.titleOfRight,(self.windowWidth-self.margin*2,150),(self.margin*2,-1),wx.ALIGN_CENTER)
		
        self.slider = wx.Slider(self.panel, (self.scaleMax-self.scaleMin)/2, (self.scaleMax-self.scaleMin)/2, self.scaleMin, self.scaleMax, pos=(self.margin, 170),
                size=(self.windowWidth-self.margin*2, -1),
                style=wx.SL_HORIZONTAL | wx.SL_AUTOTICKS   )
        self.slider.SetTickFreq((self.scaleMax-self.scaleMin)/2, (self.scaleMax-self.scaleMin)/2)


		# creates submit button which runs function saveAndClose
        if self.thisIteration==self.numOfIteration:
            self.button = wx.Button(self.panel, label="Submit",pos=(self.windowWidth/2-50,220),size=(100,50))
        else:
            self.button = wx.Button(self.panel, label="Next",pos=(self.windowWidth/2-50,220),size=(100,50))
        self.button.Bind(wx.EVT_BUTTON, self.buttonPressed)


#runs main application
        
app = wx.PySimpleApp()
frame = SliderFrame()
frame.Show()
app.MainLoop()
