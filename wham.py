import numpy as np
import glob
import csv
import math

class window:
    def __init__(self,f):
       # K=[30,50,105,131,155,171,201,221,240,258,280,309,340,365,423]
        self.center_pos = 0
        self.sc=float(f.split('_')[1].split('-')[0])
        self.framenum=f.split('_')[0].split('m')[1]
        #frm30_10-run.colvars.traj
        #frm30_10.cv
        cv_filename = 'frm{}_{}.cv'.format(self.framenum,int(self.sc))
        self.steps,self.dist=np.loadtxt(f,unpack=True)
        with open (cv_filename) as cvf:
            data = cvf.read().split('\n') #cv file split by newline
        for l in data: #per line in cv file
            words = l.split() #split each line into words
            if len(words) == 2: #if the length of words is equal to 2 (if there are two words in the line)
                if 'centers' in words[0]: #if the first word is centers
               	    self.center_pos = float(words[1]) #then the next 'word' is a float that is the center position

    def print_wind(self):
        #print('center position:', self.center_pos)
        #print('sc:', self.sc)
        #print('number of data points:', len(self.steps))
        print('f',f)
        print('steps',type(self.steps))
        print('distance',type(self.dist))
    def compute_histogram(self,bins): 
        self.occur, self.bin_edge = np.histogram(self.dist,bins)
        self.density = self.compute_density()


        #self.hist_occur, self.bin_edge = np.histogram(f,bins)



        ''' check for existence of histogram data file '''
        ''' if that is older than self.source_file, then 
        read it in and recompute hist, else just read in the hist '''

        ''' self.hist = np.hist(...) '''
    def save_histogram(self):
        file = open("histogram_{}.txt".format(self.framenum), "w")
        writer = csv.writer(file)
        
        for w in range(len(self.occur)):
            writer.writerow([self.occur[w], self.bin_edge[w]])
        file.close()
#create file names to save to and save 
        pass

#eqn 4 - for each bin compute denom by summing over all windows
    def compute_density(self): #get p_i_b
        #create new array where histogram is coverted to density by dividing each value by bin width to get density and then normal so that if you integrate you get 1 (sum bin value and bin width - same as integ). 
        results = []
        #with open("histogram_{}.txt") as csvfile:
        #reader = csv.reader(csvfile) 
        for w in range(len(self.occur)):
            #results.append([self.occur[w]/(self.bin_edge[w]-self.bin_edge[w-1])]) # divide numerator by total number of occurances in the hist
            results.append(self.occur[w]/(self.bin_edge[w]-self.bin_edge[w-1]))
            #print(results, "results_bef")
        return results
        print(results, "results_af")
        #hist_filename = "histogram_{}.txt".format(self.framenum)

  #      
   #         p_b[d]= math.exp(-beta*v[d]) #p_i_b[d]
#        self.P=

def get_bin_midpt(global_min,bin_width,bin_index):
    return global_min+(bin_index+0.5)*bin_width

def do_Wham(W,global_max,global_min,global_nb): #this function will do the wham method and W is the list of windows

    F = np.zeros(len(W)) #initialize F to 0
    
    
    nbins = global_nb
    bin_width = (global_max-global_min)/nbins 
    is_converged = False   
    while not is_converged:
        new_F = np.zeros(len(W))
        Total_P = np.zeros(nbins) #cant use len bc "object of type 'int' has no len()"
        is_converged = True
        kt=0.6
        beta = 1/kt
        denom = []
        numer = 0
        n = 0
        N = 0
        sum1 = 0
        tol = 0.000001
        for i,w in enumerate(windows):
            w.n = sum([w.occur[i] for i in range(len(w.occur))])
        #For every bin compute the denominator for each bin using equation 4
        for j in range(nbins): #b is j in eqn 4
            this_denom = 0
            #coord = W[0].bin_edge[j]-W[0].bin_edge[j-1]
            for i,w in enumerate(windows):
                v=0.5*w.sc*(get_bin_midpt(global_min,bin_width,j)-w.center_pos)**2
                #add something to catch extreme values of v to prevent floating pt error
                exp1 = math.exp(-beta*(v-F[i]))
                #print(v,"v")
                #print(F[i],"F[i]")
                if exp1 == 0: #bc initial array sets this to 0 and then gives us an error for summand (dividing by 0)
                    exp1 = 1
                #print(beta,"beta")
                #print(exp1,"exp1")
                this_denom += w.n*exp1
                #print(this_denom,"this_denom")
            denom.append(this_denom)
            #print(denom, "denom")
        for j in range(nbins-1): #if I dont put nbins-1 here, line 116 says list index for d out of range
            P_i = 0
            for w in windows:
                n = sum([w.occur[i] for i in range(len(w.occur))]) #make attribute earlier in hist
                #print(w.density[j], "w.density[j]")
                #print(type(w.density[0]), "type(w.density[j])")
                #print("range",range(nbins))
                d = w.density[j]

                summand = n*d/denom[j]
                #print(denom[j], "denom[j]")
                #print(summand, "summand")
                P_i += summand
                #print(P_i, "P_i")
            Total_P[j] = P_i #eqn 4
   #eqn 5 below 
           
        for i,w in enumerate(windows):
            for j in range(nbins):
                v=0.5*w.sc*(get_bin_midpt(global_min,bin_width,j)-w.center_pos)**2 #bias potential
                sum1 += Total_P[j]*math.exp(-beta*v) # v is dependent on b
            F_i = (-1/beta)*math.log(sum1)
            new_F[i] = F_i
        #do a sum of sq differences element by element
        ssqd = np.sum((F-new_F)**2)
        print(ssqd, "ssqd")
        if ssqd>tol:
            is_converged = False
            F = new_F
    return Total_P

if __name__ == '__main__':

    global_min = 40.0  #?
    global_max = 60.0  # ?
    global_nb = 100
    
    ''' expects *traj files in directory are complete set '''
    filenames=glob.glob('*.traj') #instead of glob loop over k in traj name

    windows = [] #empty list called 'windows'
    for f in filenames:
        this_window = window(f) #here we are calling the class window as 'this_window' variable for file f
        #this_window.print_wind() #here we are calling the def print_wind() function in the class
        windows.append(this_window) #add
        #print(f,'has {:d} lines of data for spring constant {:.2f}'.format(len(steps),k))
    bins = np.linspace(global_min, global_max, global_nb)
    for w in windows:
        w.compute_histogram(bins)
        w.save_histogram()
    do_Wham(windows,global_max,global_min,global_nb)
    
  #  P_unbiased = do_Wham(windows)

