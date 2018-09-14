import numpy as np
import trident
import yt
import os
import sys
import matplotlib.pyplot as plt
from quasar_scan import *
from parse_vela_metadata import Rdict, Ldict

#precondition: assumes there are only two levels of depth within the output folder
#postcondition: returns a list of all textfiles
def get_all_textfiles(inquasarscan = True):
    
    #pathname by default starts in output
    path = "output2.0"
    if not inquasarscan:
        path = "quasarscan/output2.0"
    textfiles = []
    
    #gets all folders in output
    dirs = os.listdir(path)
    for folderName in dirs:
        if not folderName.startswith("."):
            folderPath = path + "/" + folderName
            folderDirs = os.listdir(folderPath)
            for fileName in folderDirs:
                if not fileName.startswith("."):
                    textfiles.append(os.path.join(folderPath,fileName))
    print ("All textfiles loaded!")
    return textfiles

def get_VELA_folder_textfiles(folderName):
    
    #pathname by default starts in output
    path = "output2.0/" + folderName + "coldensinfo"
    
    textfiles = []
    
    #gets all folders in output
    dirs = os.listdir(path)
    for fileName in dirs:
        folderPath = path + "/"
        textfiles.append(os.path.join(path,fileName))
    return textfiles
    
def sort_ions(ions,flat = True):
    def sort_ions_one_element(ions,element):
        from trident import from_roman,to_roman
        nums = [None]*len(ions)
        toreturn = []
        for i in range(len(ions)):
            nums[i] = from_roman(ions[i].split(" ")[1])
        nums.sort()
        for val in nums:
            toreturn.append("%s %s"%(element,to_roman(val)))
        return toreturn
    ions = list(ions)
    ions.sort()
    index = 0
    element = ions[index].split(" ")[0]
    tosort = []
    toreturn = []
    while index < len(ions):
        if ions[index].split(" ")[0] == element:
            tosort.append(ions[index])
            index += 1
        else:
            toreturn.append(sort_ions_one_element(tosort,element))
            element = ions[index].split(" ")[0]
            tosort = []
    toreturn.append(sort_ions_one_element(tosort,element))
    if flat:
        toreturn = [item for sublist in toreturn for item in sublist]
    return toreturn

stringcriteria = ["ions","version","simname","simnum","has_intensives"]
intensives = ["Z","T","n"]
intensiveslabels = {"Z":"avg log metallicity","T":"avg log temperature","n":"avg log density"}
intensivespositions = {"Z":-1,"T":-2,"n":-3}
sightline_xVars = ["r","rdivR","theta","phi"]
param_xVars = ["redshift","a0","Mvir","gas_Rvir","star_Rvir","dm_Rvir","sfr","ssfr","L_mag"]

class MultiQuasarSpherePlotter():
    #USER MUST EXPLICITLY CALL GET_QUASAR TO INPUT INTO ALL OTHER MULTIQUASARSPHEREPLOTTER METHODS AS NECESSARY
    
    
    #param: textfiles     if a list of textfiles is specified, those specific textfiles will be loaded; else,
    #                     all textfiles in output are loaded
    def __init__(self, textfiles = None, cleanup = False,plots = "mean"):
        self.plots = "mean"
        self.avgfn = np.mean
        self.setPlots(plots)
        self.quasarLineup = []
        if textfiles is None:
            textfiles = get_all_textfiles()
        for textfile in textfiles:
            try:
                simparams,scanparams,ions,data,gasbins = read_values(textfile)
                q = QuasarSphere(simparams = simparams,scanparams = scanparams,ions = ions,data = data, gasbins = gasbins, readonly = True)
                if self.pass_safety_check(q):
                    self.quasarLineup.append(q)
                elif cleanup:
                    todo = raw_input("file %s did not pass safety check. Remove it? (y/n)"%textfile).lower()
                    os.remove(textfile) if todo == 'y' else None
            except Exception as e:
                print(textfile + " could not load because:")
                print(e)
        self.quasarArray = np.array(self.quasarLineup)
        self.currentQuasarArray = []
        for q in self.quasarArray:
            self.currentQuasarArray.append(q)
        self.currentQuasarArray = np.array(self.currentQuasarArray)
        self.currentQuasarArrayName = ''
        
        if len(self.currentQuasarArray) == 0:
            print ("There are no quasarspheres stored in currentQuasarArray!")

    #tests each condition, each condition is a safety check
    def pass_safety_check(self, q):
        minlength = 1
        minions = []
        if q.length_reached < minlength:
            print "Length for %s is not valid." %(q.simname + "z" + str(q.rounded_redshift))
            return False
        elif not all(x in q.ions for x in minions):
            print "Not all necessary ions present in %s."%(q.simname + "z" + str(q.rounded_redshift))
            return False
        try:
            q.final_a0
        except: 
            print "metadata for %s is not valid." %(q.simname + "z" + str(q.rounded_redshift))
            return False
        return True
        
    def reset_current_Quasar_Array(self):
        self.currentQuasarArray = []
        for q in self.quasarArray:
            self.currentQuasarArray.append(q)
        self.currentQuasarArray = np.array(self.currentQuasarArray)
        self.currentQuasarArrayName = ''
    
    #param bins either serves as an array with each element being as a cutpoint, or a single value
    #follows array indexing exclusivity and inclusivity rules
    def sort_by(self, criteria, bins = [0,np.inf], reset = False, exploration_mode = False,atEnd = False,onlyNonempty = False,splitEven = 0):
        if not (criteria in self.currentQuasarArray[0].__dict__.keys()):
            print ("Criteria " + criteria + " does not exist. Please re-enter a valid criteria.")
            return
        elif criteria == "ions":
            print("You cannot sort by 'ions'")
            return
        if criteria in self.currentQuasarArrayName:
            print("Already constrained by %s. Please reset instead of further constraining."%criteria)
            return
        if isinstance(bins, float) or isinstance(bins, int) or \
            (len(bins) == 1 and isinstance(bins[0], float)) or (len(bins) == 1 and isinstance(bins[0], int)):
            if isinstance(bins, list) or isinstance(bins,np.ndarray):
                bins = bins[0]
            bins = np.array([0.0, bins, np.inf])
        elif isinstance(bins, str) and criteria in stringcriteria:
             bins = [bins]
        sorter = MultiSphereSorter(self.currentQuasarArray,exploration_mode = exploration_mode)
        if splitEven:
            labels, bins, quasarBins = sorter.splitEven(criteria,splitEven,atEnd = atEnd)
        else:
            labels, bins, quasarBins = sorter.sort(criteria,bins,atEnd = atEnd)
        if quasarBins is None:
            return
        if reset == True:
            self.reset_current_Quasar_Array()
        empty = True
        nonemptyArray = []
        nonemptyLabelArray = []
        for i in range(len(quasarBins)):
            item = quasarBins[i]
            if len(item)>0:
                nonemptyArray.append(item)
                nonemptyLabelArray.append(labels[i])
                empty = False
        if onlyNonempty:
            return np.array(nonemptyLabelArray),bins, np.array(nonemptyArray)
        print "Bins are empty." if empty else ""
        return labels,bins, quasarBins
        
    def constrain_current_Quasar_Array(self, constrainCriteria, bins, exploration_mode = False,atEnd = False,splitEven = None,extra_title = ""):
        if len(self.currentQuasarArray)>0 and not (constrainCriteria in self.currentQuasarArray[0].__dict__.keys()):
            print ("Constrain criteria " + constrainCriteria + " does not exist. Please re-enter a valid criteria.")
            return
        if constrainCriteria in self.currentQuasarArrayName:
            print("Already constrained by %s. Please reset instead of further constraining."%constrainCriteria)
            return
        if isinstance(bins, list):                
            if len(bins) != 2 and constrainCriteria not in stringcriteria:
                print ("Length of bins must be 2: [lower,upper]")
                return
        elif isinstance(bins,str) and constrainCriteria in stringcriteria:
            bins = [bins]
        sorter = MultiSphereSorter(self.currentQuasarArray,exploration_mode = exploration_mode)
        if splitEven is None:
            labels, bins, temp = sorter.sort(constrainCriteria,bins,atEnd = atEnd)
        else:
            labels,bins,temp = sorter.splitEven(constrainCriteria,2,atEnd = atEnd)
            if splitEven == "high":
                take = 1
            elif splitEven == "low":
                take = 0
            else:
                print("please use splitEven = 'low' or 'high'")
                return
            labels = np.array([labels[take]])
            bins = np.array([bins[take],bins[take+1]])
            temp = np.array([temp[take]])
        if temp is None:
            return
        
        self.currentQuasarArray = np.unique(np.concatenate(temp))
        
        self.currentQuasarArrayName += constrainCriteria 
        if not constrainCriteria in stringcriteria:
            if bins[0] == 0.0:
                lowlabel = "lessthan"
            elif bins[0] > 100 or bins[0]<.1:
                lowlabel = "%1.1e"%bins[0]
            else:
                lowlabel = "%1.1f"%bins[0]
            if bins[1] == np.inf:
                highlabel = "andhigher"
            elif bins[1] > 100 or bins[1]<.1:
                highlabel = "-%1.1e"%bins[1]
            else:
                highlabel = "-%1.1f"%bins[1]
            self.currentQuasarArrayName += "%s%s"%(lowlabel,highlabel)
        else:
            for acceptedValue in bins:
                self.currentQuasarArrayName += acceptedValue.replace(" ","")
        return bins
    #summary: plots an a pyplot errorbar graph with the x-axis being either theta, phi, or the radius; 
    #         the y-axis points are the mean column densities with a spread of +/- the error
    #         quasarArray is the array of quasars to be plotted


    def setPlots(self,plots):
        def getquartiles(data):
            return np.array([np.median(data) - np.percentile(data,40),np.percentile(data,60)-np.median(data)])
        def getstderr(data):
            return np.array([np.std(data)/np.sqrt(len(data)),np.std(data)/np.sqrt(len(data))])
        if plots == "mean":
            self.plots = "mean"
            self.avgfn = np.mean
            self.errfn = getstderr
        elif plots == "median":
            self.plots = "median"
            self.avgfn = np.median
            self.errfn = getquartiles
        elif plots == "med_noquartiles":
            self.plots = "median_std"
            self.avgfn = np.median
            self.errfn = getstderr

    def get_xy_type0(self,xVar,ions,quasarArray,rlims):
        quasarArray = [self.currentQuasarArray]
        xs = np.empty(len(ions),dtype = object)
        ys = np.empty(len(ions),dtype = object)
        for i in range(len(ions)):
            ion = ions[i]
            if xVar in sightline_xVars:
                x,y = self.get_xy_type1(xVar,ion,quasarArray,rlims)
            elif xVar in param_xVars:
                x,y = self.get_xy_type2(xVar,ion,quasarArray,rlims)
            xs[i] = x[0]
            ys[i] = y[0]
        return xs,ys

    def get_xy_type1(self,xVar,ion,quasarArray,rlims):
        if rlims is None:
            rlims = [0.1,np.inf]
        vardict = {"theta":1,"phi":2,"r":3,"rdivR":3}
        distances = "kpc" if xVar == "r" else "Rvir"
        gqary = []
        xs = np.empty(len(quasarArray),dtype = object)
        ys = np.empty(len(quasarArray),dtype = object)
        for i in range(len(quasarArray)):
            gq = GeneralizedQuasarSphere(quasarArray[i],distance=distances)
            if gq.number == 0:
                xs[i] = np.empty(0)
                ys[i] = np.empty(0)
                continue
            xs[i] = gq.info[:,vardict[xVar]]
            ys[i] = gq.info[:,gq.get_ion_column_num(ion)]
            rs = gq.info[:,3]
            acceptedLines = np.logical_and(rlims[0]<=rs,rs<=rlims[1])
            xs[i] = xs[i][acceptedLines]
            ys[i] = ys[i][acceptedLines]
            xs[i] = xs[i][ys[i]>0]
            ys[i] = ys[i][ys[i]>0]
            ys[i] = np.log10(ys[i])
        return xs,ys
    
    def get_xy_type2(self,xVar,ion,quasarArray,rlims):
        if rlims is None:
            rlims = np.array([0.1,1.0])
        xs = np.empty(len(quasarArray),dtype = object)
        ys = np.empty(len(quasarArray),dtype = object)
        for i in range(len(quasarArray)):
            ary = quasarArray[i]
            xs[i] = np.empty(0)
            ys[i] = np.empty(0)
            for q in ary:
                x = eval("q."+xVar)
                cdens = q.info[:,q.get_ion_column_num(ion)]
                rs = q.info[:,3]
                acceptedLines = np.logical_and(rlims[0]*q.Rvir<=rs*q.code_unit_in_kpc,\
                                               rs*q.code_unit_in_kpc<=rlims[1]*q.Rvir)
                cdens = cdens[acceptedLines]
                cdens = cdens[cdens>0]
                xs[i] = np.append(xs[i],np.tile(x,len(cdens)))
                ys[i] = np.append(ys[i],np.log10(cdens))
        return xs,ys
    
    def get_xy_type3(self,xVar,yVar,quasarArray):
        xs = np.empty(len(quasarArray),dtype = object)
        ys = np.empty(len(quasarArray),dtype = object)
        for i in range(len(quasarArray)):
            ary = quasarArray[i]
            xs[i] = np.empty(0)
            ys[i] = np.empty(0)
            for q in ary:
                x = eval("q."+xVar)
                y = eval("q."+yVar)
                xs[i] = np.append(xs[i],x)
                ys[i] = np.append(ys[i],y)
        return xs,ys

    def combine_xs(self,x_variable,tolerance):
        x_values_not_averaged = np.unique(x_variable)
        x_values = []
        i = 0
        while i < len(x_values_not_averaged)-1:
            currentList = [x_values_not_averaged[i]]
            numtocombine=1
            while i+numtocombine < len(x_values_not_averaged) and \
                    x_values_not_averaged[i+numtocombine]/x_values_not_averaged[i]-1 <= tolerance:
                currentList.append(x_values_not_averaged[i+numtocombine])
                numtocombine+=1
            i+=numtocombine
            x_values.append(self.avgfn([currentList]))
        if i == len(x_values_not_averaged)-1:
            x_values.append(x_values_not_averaged[-1]) 
        return np.array(x_values)

    def get_savefig_name(self,ion,labels,xVar,plot_type):
        if plot_type in [0]:
            ionNameNoSpaces = str(labels).strip("[]").replace("'","").replace(" ","")
            name = "%s_ErrorBar_%s_%s_%s" % (self.currentQuasarArrayName, ionNameNoSpaces, xVar, self.plots)
        elif plot_type in [1,2,3]:
            ionNameNoSpaces = ion.replace(" ","")
            binVariables = labels[0].split(" ")
            for binVariable in binVariables:
                if binVariable in ["<",">"]:
                    continue
                try:
                    number = float(binVariable)
                    continue
                except:
                    break
            name = "%s_ErrorBar_%s_%s_%s_%s" % (self.currentQuasarArrayName, ionNameNoSpaces, xVar, binVariable, self.plots)
        if name[0] == "_":
            name = name[1:]
        return name

    def plot_err(self, ion, quasarArray = None, xVar = "r", save_fig = False, \
                 reset = False, labels = None,extra_title = "",rlims = None,\
                 tolerance = 1e-5,dots = False,logx = False,average = None):
        print("Current constraints (name): "+self.currentQuasarArrayName)
        plt.figure()
        if isinstance(ion,list):
            ions = ion
            plot_type = 0
            for ion in ions:
                assert not ion in intensives
            assert quasarArray is None
            assert labels is None
            labels = ions
            xarys,yarys = self.get_xy_type0(xVar,ions,quasarArray,rlims)
        elif not quasarArray is None and xVar in sightline_xVars:
            plot_type = 1
            assert isinstance(ion,str)
            assert len(labels) == len(quasarArray)
            xarys,yarys = self.get_xy_type1(xVar,ion,quasarArray,rlims)
        elif (not quasarArray is None) and xVar in param_xVars and (not ion in param_xVars):
            plot_type = 2
            assert isinstance(ion,str)
            assert len(labels) == len(quasarArray)
            xarys,yarys = self.get_xy_type2(xVar,ion,quasarArray,rlims)
        elif ion in param_xVars and xVar in param_xVars:
            yVar = ion
            plot_type = 3
            assert not quasarArray is None
            assert len(labels) == len(quasarArray)
            if not dots is None: 
                dots = True
            xarys,yarys = self.get_xy_type3(xVar,yVar,quasarArray)
        else:
            print("Requirements not met (type could not be detected). Cannot plot.")
            return

        if not average is None:
            self.setPlots(average)
        xs = np.empty(len(yarys),dtype = object)
        ys = np.empty(len(yarys),dtype = object)
        yerrs = np.empty(len(yarys),dtype = object)
        for i in range(len(yarys)):
            unique_xs = self.combine_xs(xarys[i],tolerance)
            xs[i] = unique_xs
            ys[i] = np.zeros(len(unique_xs))
            yerrs[i] = np.zeros((len(unique_xs),2))
            mask = np.ones(len(unique_xs),dtype = bool)
            for j in range(len(unique_xs)):
                x_value = unique_xs[j]
                yvals = yarys[i][abs(xarys[i]/x_value - 1)<=tolerance]
                if len(yvals) == 0:
                    mask[j] = False
                    continue
                ys[i][j] = self.avgfn(yvals)
                yerrs[i][j] = self.errfn(yvals)
            xs[i] = xs[i][mask]
            ys[i] = ys[i][mask]
            yerrs[i] = yerrs[i][mask]
            
        sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
                   "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)","phi" \
                   :"azimuthal viewing angle (rad)"}
        param_unit_labels = {"redshift":"z","a0":"a","Mvir":"Virial Mass (Msun)",\
                            "gas_Rvir":"Gas Mass within Rvir (Msun)","star_Rvir":"Stellar Mass within Rvir (Msun)",\
                            "dm_Rvir":"Dark Matter Mass within Rvir (Msun)","sfr":"Star Formation Rate (Msun yr-1)",\
                            "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)","L_mag":"Magnitude of Angular Momentum"}

        if xVar in sightline_xVars:
            plt.xlabel(sightline_unit_labels[xVar])
        elif xVar in param_xVars:
            plt.xlabel(param_unit_labels[xVar])
        if ion in intensives:
            plt.ylabel(intensiveslabels[ion])
            cd = ""
        elif ion in param_xVars:
            plt.ylabel(param_unit_labels[yVar])
            cd = ""
        else:
            plt.ylabel("log col dens")
            cd = "Column Density "
        if plot_type in [0]:
            ionstr = str(ions).strip("[]").replace("'","")
        elif plot_type in [1,2,3]:
            ionstr = ion
        plt.title('%s %sAverages (%s) %s'%(ionstr, cd, self.plots, extra_title))

        for i in range(len(xs)):
            x = xs[i]
            y = ys[i]
            yerr = np.transpose(yerrs[i])
            label = labels[i] 
            #change fmt to . or _ for a dot or a horizontal line
            fmtdict = {"mean":',',"median_std":',',"median":"."}
            if dots:
                plt.plot(x,y,"o",label = label)
            else:
                plt.errorbar(x,y, yerr=yerr, fmt=fmtdict[self.plots], capsize = 3, label = label)
        plt.legend()
        
        if logx:
            plt.xscale('log')

        if reset:
            self.reset_current_Quasar_Array()
            
        if save_fig:
            name = self.get_savefig_name(ion,labels,xVar,plot_type)
            plt.savefig("plots/"+name + ".png")
            return plt,"plots/"+name + ".png"
        return plt

    
    #summary: plots histogram(s) of ion column density vs incrementing xvariable, 
    #         uses a color bar to show the percentage of sightlines for a certain column density at a specific x value
    def plot_hist(self, q, simname = None,xvariable = "r",zeros = "ignore",\
                  weights = True,save_fig = None,ns = (42,15),do_ions = "all"):
        if not simname:
            simname = q.simparams[0]
        if do_ions == "all":
            ions = q.ions
        else:
            ions = do_ions
   
        if xvariable in ["r","r>0"]:
            conversion = q.simparams[10]
        elif xvariable in ["rdivR","rdivR>0"]:
            if q.simparams[5] > 0:
                conversion = q.simparams[10]/q.simparams[5]
            else:
                print("No virial radius found")
                return
        else:
            if q.simparams[5] < 0:
                print("No metadata, angle plots will be arbitrary axis.")
            conversion = 1

        vardict = {"theta":1,"phi":2,"r":3,"r>0":3,"rdivR":3,"rdivR>0":3}
        if not (xvariable in vardict):
            print ("Inputted xvariable not found. Please enter a valid xvariable.")
            return

        xVarsArray = q.info[:, vardict[xvariable]]*conversion
        #ion,xvars,cdens,simname
        for i in range(len(q.ions)):
            end = q.scanparams[6]
            if q.ions[i] in ions:
                plot2dhist(q.ions[i],xVarsArray,\
                       q.info[:end,11+i],simname,xvariable = xvariable, ns = ns,zeros = zeros,\
                       weights = weights,save_fig = save_fig,z = q.simparams[1])

#summary: helper method for plot_hist                
def plot2dhist(ion,xvars,cdens,simname,xvariable = "r",ns = (42,15),zeros = "ignore",weights = True, save_fig = None, z = None):
    if zeros == "ignore":
        xvars = xvars[cdens>0]
        cdens = cdens[cdens>0]
        logdens = np.log10(cdens)
    else:
        logdens = np.log10(np.maximum(cdens,1e-15))
    if xvariable == "r>0" or xvariable == "rdivR>0":
        logdens = logdens[xvars>0.0]
        xvars = xvars[xvars>0.0]
    nx = ns[0]
    ny = ns[1]
    plotvars = {"r":"r","r>0":"r","rdivR":"r","rdivR>0":"r","theta":"theta","phi":"phi"}
    plotvar = plotvars[xvariable]
    if weights:
        weight = xvars*0.0
        for i in range(len(xvars)):
            weight[i] = 1.0/len(xvars[xvars==xvars[i]])
        H, xedges, yedges = np.histogram2d(xvars, logdens, bins=[nx,ny],weights = weight)
        cbarlabel = "Fraction of lines for fixed %s"%(plotvar)
    else:
        H, xedges, yedges = np.histogram2d(xvars, logdens, bins=[nx,ny])
        cbarlabel = "Total number of lines"
    H = H.T
    X, Y = np.meshgrid(xedges, yedges)
    plt.pcolormesh(X,Y, H, cmap=hotcustom)
    plt.title("distribution of "+ion+" in "+simname+" at z="+str(z)[:4])
    # set the limits of the plot to the limits of the data
    #plt.axis([x.min(), x.max(), y.min(), y.max()])
    plt.colorbar(label = cbarlabel)
    x1,x2,y1,y2 = plt.axis()
    dx = x2-x1
    dy = y2-y1
    plt.axis((x1-dx*0.1,x2+dx*0.1,y1-dy*0.1,y2+dy*0.1))
    xlabels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir","theta":"viewing angle (rad)","phi":"azimuthal viewing angle (rad)"}
    plt.xlabel(xlabels[xvariable])
    plt.ylabel("log col dens")
    if save_fig:
        if save_fig == "default" :
            save_fig = simname + "_" + plotvar + "_z" +str(z)[:4]
        name = save_fig+"_"+ion.replace(" ","")
        if weights:
            name +="_w"
        if zeros == "ignore":
            name +="_nozeros"
        plt.savefig(name+".png")
    plt.show()

class MultiSphereSorter(object):
    def __init__(self,myArray,exploration_mode = False):
        self.array = myArray
        self.exploration_mode = exploration_mode
    #param bins the lower parameter is inclusive, the upper parameter is exclusive
    def sort(self, criteria, bins,atEnd = False):
        labels = self.make_labels(criteria, bins,atEnd = atEnd)
        if criteria in stringcriteria:
            sortfn = self.sort_by_strparam
        else:
            sortfn = self.sort_by_default
        while self.exploration_mode:
            fakeBins = sortfn(criteria, bins, atEnd = atEnd)
            print "Bins will be categorized in the following structure: \n"
            for index in range(len(fakeBins)):
                oneBin = fakeBins[index]
                print (labels[index] + " has " + str(len(oneBin)) + " elements out of " + str(len(self.array)) + "\n")
            response = raw_input("Continue? ([Y]/N) or enter new 'bin' parameter.\n")
            if response.lower() == "n":
                return None, None, None
            elif response.lower() == "y" or response == "":
                self.exploration_mode = False
            else:
                while True:
                    try:
                        bins = eval(response)
                        if isinstance(bins, float) or isinstance(bins, int) or \
                            (len(bins) == 1 and isinstance(bins[0], float)) or (len(bins) == 1 and isinstance(bins[0], int)):
                            if isinstance(bins, list) or isinstance(bins,np.ndarray):
                                bins = bins[0]
                            bins = np.array([0.0, bins, np.inf])
                        elif isinstance(bins, str) and criteria in stringcriteria:
                            bins = [bins]
                        break
                    except Exception as e:
                        print e
                        response = raw_input("Could not evaluate %s. Please re-enter.\n"%response)       
            labels = self.make_labels(criteria, bins, atEnd = atEnd)
        return labels, bins, sortfn(criteria, bins, atEnd = atEnd)

    def sort_by_strparam(self, criteria, acceptedValues,atEnd = False):
        if atEnd:
            print("dont use atEnd for string sort")
            return
        criteriaArray = self.get_criteria_array(criteria)
        resArray = np.empty(len(acceptedValues),dtype = 'object')
        for i in range(len(acceptedValues)):
            toAdd = []
            for j in range(len(criteriaArray)):
                add = False
                if criteria == "ions" and acceptedValues[i] in criteriaArray[j]:
                    add = True
                elif acceptedValues[i] == criteriaArray[j]:
                    add = True
                if add:
                    toAdd.append(self.array[j])
            resArray[i] = np.array(toAdd)
        return np.array(resArray)   
    
    '''These are all the same as sort_by_default
    def sort_by_redshift(self, mq, criteria, bins):
    def sort_by_rmax(self, mq, criteria, bins):
    def sort_by_a0(self, mq, criteria, bins):
    def sort_by_length(self, mq, criteria, bins):
    def sort_by_R(self, mq, criteria, bins):    
    def sort_by_len_th_arr(self, mq, criteria, bins):    
    def sort_by_len_r_arr(self, mq, criteria, bins):
    def sort_by_len_phi_arr(self, mq, criteria, bins):
    def sort_by_Rvir(self, mq, criteria, bins):
    def sort_by_gas_Rvir(self, mq, criteria, bins):
    def sort_by_dm_Rvir(self, mq, criteria, bins):
    def sort_by_star_Rvir(self, mq, criteria, bins):
    def sort_by_sfr(self, mq, criteria, bins):
    '''
    
    def sort_by_default(self, criteria, bins, atEnd = False):
        if len(bins) == 1:
            criteriaArray = self.get_criteria_array(criteria, atEnd = atEnd)
            resArray = []
            for index in range(len(criteriaArray)):
                intersection = np.intersect1d(bins, criteriaArray[index])
                if len(intersection) > 0:
                    resArray.append(self.array[index]) 
            return resArray 
        resultBins = [None] * (len(bins)-1)
        criteriaArray = self.get_criteria_array(criteria,atEnd=atEnd)
        for index in range(len(bins)-1):
            booleanindices = np.logical_and(criteriaArray >= bins[index], criteriaArray < bins[index+1]) 
            toadd = self.array[booleanindices]
            resultBins[index] = toadd
        return np.array(resultBins)
    
    def splitEven(self,criteria,num,atEnd = False):
        if criteria in stringcriteria:
            print("cannot splitEven over string criteria")
        criteriaArray = self.get_criteria_array(criteria, atEnd = atEnd)
        sortedcriteriaArray = np.sort(criteriaArray)
        quotient = len(sortedcriteriaArray) // num
        if quotient == 0:
            print("Warning: Number of bins exceeds length of criteria array. Not all bins will be filled.")
        remainder = len(sortedcriteriaArray) % num
        bin_edges = np.zeros(num + 1)
        bin_edges[0] = 0
        bin_edges[-1] = np.inf
        j = 0
        for i in range(1,num):
            if i <= remainder:
                bin_edges[i] = np.mean(sortedcriteriaArray[i*quotient+j : i*quotient+2+j])
                j+=1
            else:
                bin_edges[i] = np.mean(sortedcriteriaArray[i*quotient-1+j : i*quotient+1+j])
        return self.sort(criteria, bin_edges,atEnd = atEnd)
    
    
    def get_criteria_array(self, criteria,atEnd = False):
        if atEnd == False:
            res = []
            for q in self.array:
                criteria_vals = eval("q." + criteria)
                res.append(criteria_vals)
            return np.array(res)
        else:
            quasar_array = self.array
            final_a = np.zeros(len(quasar_array))
            for i in range(len(quasar_array)):
                final_a[i] = quasar_array[i].final_a0
            min_a = "%1.3f"%min(final_a)
            criteria_final = np.zeros(len(quasar_array))
            for index in range(len(quasar_array)):
                q = quasar_array[index]
                criteria_final[index] = q.get_criteria_at_a(min_a,criteria)
            return criteria_final
    
    def make_labels(self,criteria, bins, atEnd=False):
        labels = []
        if criteria in stringcriteria:
            labels = bins
        else:
            if atEnd:
                quasar_array = self.array
                final_a = np.zeros(len(quasar_array))
                for i in range(len(quasar_array)):
                    final_a[i] = quasar_array[i].final_a0
                min_a = "%1.3f"%min(final_a)
                criteria = criteria + "_"+min_a 
            for index in range(len(bins)-1):
                low = bins[index]
                high = bins[index+1]
                lowstr = "%.1e"%low if low < 0.1 or low > 100.0 else str(low)[:4]
                highstr = "%.1e"%high if high < 0.1 or high > 100.0 else str(high)[:4]
                if low == 0.0:
                    uniqueName = "%s < %s"%(criteria, highstr)
                elif high == np.inf:
                    uniqueName = "%s > %s"%(criteria, lowstr)
                else:
                    uniqueName = "%s < %s < %s"%(lowstr,criteria,highstr)
                labels.append(uniqueName)
        return labels
