import numpy as np
import trident
import yt
import os
import sys
import matplotlib.pyplot as plt
from quasar_scan import *
import gasbinning
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
sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
           "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)","phi" \
           :"azimuthal viewing angle (rad)"}
param_unit_labels = {"redshift":"z","a0":"a","Mvir":"Virial Mass (Msun)",\
                    "gas_Rvir":"Gas Mass within Rvir (Msun)","star_Rvir":"Stellar Mass within Rvir (Msun)",\
                    "dm_Rvir":"Dark Matter Mass within Rvir (Msun)","sfr":"Star Formation Rate (Msun yr-1)",\
                    "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)","L_mag":"Magnitude of Angular Momentum"}

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
    def sort_by(self, criteria, bins = [0,np.inf], reset = False, exploration_mode = False,\
        atEnd = False,onlyNonempty = False,splitEven = 0,reverse = False):
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
        if reverse:
            def reversearray(ary):
                ary = list(ary)
                ary.reverse()
                return np.array(ary)
            labels = reversearray(labels)
            bins = reversearray(bins)
            quasarBins = reversearray(quasarBins)
        return labels,bins, quasarBins
    
    def constrain_via_gasbins(self,gasbintype=None):
        if gasbintype == None:
            gasbintype = raw_input("Available bins are: %s"%gasbinning.possible_bin_types)
        g = gasbinning.GasBinsHolder(bins=[gasbintype])
        toReturn = []
        for q in self.currentQuasarArray:
            if g.get_bin_str() in q.gasbins.get_bin_str():
                toReturn.append(q)
        self.currentQuasarArray = toReturn
    
    def constrain_current_Quasar_Array(self, constrainCriteria, bins, exploration_mode = False,atEnd = False,splitEven = None):
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


    def setPlots(self,plots,quartiles = None):
        if isinstance(plots,tuple):
            if len(plots) == 3:
                quartiles = (plots[1],plots[2])
            elif len(plots) == 2:
                if isinstance(plots[1],tuple):
                    quartiles = plots[1]
                else:
                    quartiles = (50-plots[1],50+plots[1])
            plots = plots[0]
        if not quartiles:
            quartiles = (40,60)
        def getquartiles(data):
            return np.array([np.nanmedian(data) - np.nanpercentile(data,quartiles[0]),\
                             np.nanpercentile(data,quartiles[1])-np.nanmedian(data)])
        def getstderr(data):
            l = len(data[~np.isnan(data)])
            return np.array([np.nanstd(data)/np.sqrt(l),np.nanstd(data)/np.sqrt(l)])
        if plots == "mean" or plots == "std":
            self.plots = "mean"
            self.avgfn = np.nanmean
            self.errfn = getstderr
        elif plots == "median" or plots == "med":
            self.plots = "median"
            self.avgfn = np.nanmedian
            self.errfn = getquartiles
        elif plots == "med_noquartiles" or plots == "median_std":
            self.plots = "median_std"
            self.avgfn = np.nanmedian
            self.errfn = getstderr
        elif plots == "scatter":
            self.plots = "scatter"
            self.avgfn = np.nanmean
            self.errfn = getstderr
    
    def get_yVar_from_str(self,gq,stringVar):
        def split_by_ops(s):
            s = s.replace("(","_splitchar_")
            s = s.replace(")","_splitchar_")
            s = s.replace("/","_splitchar_")
            s = s.replace("*","_splitchar_")
            s = s.replace("+","_splitchar_")
            #s = s.replace("-","_splitchar_")
            return filter(None,s.split("_splitchar_"))
        def bylength(word1,word2):
            return len(word2)-len(word1)
        def sortlist(a):
            a.sort(cmp=bylength)
            return a
        strings_to_find = split_by_ops(stringVar)
        new_str_to_eval = stringVar
        strings_to_replace_with = {}
        for i,s in enumerate(strings_to_find):
            try:
                eval(s)
                string_to_replace_with = s
            except:
                string_to_replace_with = "gq.info[:,%d]"%gq.get_ion_column_num(s)
            strings_to_replace_with[s] = string_to_replace_with
        for s in sortlist(strings_to_replace_with.keys()):
            new_str_to_eval = new_str_to_eval.replace(s,strings_to_replace_with[s])
        to_return = eval(new_str_to_eval)
        return to_return

    def get_xy_type0(self,xVar,yVars,quasarArray,rlims):
        quasarArray = [self.currentQuasarArray]
        xs = np.empty(len(yVars),dtype = object)
        ys = np.empty(len(yVars),dtype = object)
        for i in range(len(yVars)):
            yVar = yVars[i]
            if xVar in sightline_xVars:
                x,y = self.get_xy_type1(xVar,yVar,quasarArray,rlims)
            elif xVar in param_xVars:
                x,y = self.get_xy_type2(xVar,yVar,quasarArray,rlims)
            xs[i] = x[0]
            ys[i] = y[0]
        return xs,ys

    def get_xy_type1(self,xVar,yVar,quasarArray,rlims):
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
            ys[i] = self.get_yVar_from_str(gq,yVar)
            rs = gq.info[:,3]
            acceptedLines = np.logical_and(rlims[0]<=rs,rs<=rlims[1])
            xs[i] = xs[i][acceptedLines]
            ys[i] = ys[i][acceptedLines]
        return xs,ys
    
    def get_xy_type2(self,xVar,yVar,quasarArray,rlims):
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
                y = self.get_yVar_from_str(q,yVar)
                rs = q.info[:,3]
                acceptedLines = np.logical_and(rlims[0]*q.Rvir<=rs*q.code_unit_in_kpc,\
                                               rs*q.code_unit_in_kpc<=rlims[1]*q.Rvir)
                y = y[acceptedLines]
                xs[i] = np.append(xs[i],np.tile(x,len(y)))
                ys[i] = np.append(ys[i],y)
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
    
    def get_xy_type4(self, xVar, yVar, quasarArray, rlims):
        if rlims is None:
            rlims = [0.1,1.0]
        xs = np.empty(len(quasarArray),dtype = object)
        ys = np.empty(len(quasarArray),dtype = object)
        for i in range(len(quasarArray)):
            xs[i] = np.empty(len(quasarArray[i]),dtype = object)
            ys[i] = np.empty(len(quasarArray[i]),dtype = object)
            for j in range(len(quasarArray[i])):
                q = quasarArray[i][j]
                x = self.get_yVar_from_str(q,xVar)
                y = self.get_yVar_from_str(q,yVar)
                rs = q.info[:,3]
                acceptedLines = np.logical_and(rlims[0]*q.Rvir<=rs*q.code_unit_in_kpc,\
                                               rs*q.code_unit_in_kpc<=rlims[1]*q.Rvir)
                x = x[acceptedLines]
                y = y[acceptedLines]
                xs[i][j] = x
                ys[i][j] = y
        return xs,ys

    def values_within_tolerance(self,x_in_list,x_comp,tolerance):
        if x_comp == 0:
            return np.abs(x_in_list)<=tolerance
        else:
            return x_in_list/x_comp-1<=tolerance

    def combine_xs(self,x_variable,tolerance):
        x_values_not_averaged = np.unique(x_variable)
        x_values = []
        i = 0
        while i < len(x_values_not_averaged)-1:
            currentList = [x_values_not_averaged[i]]
            numtocombine=1
            while i+numtocombine < len(x_values_not_averaged) and \
                    self.values_within_tolerance(x_values_not_averaged[i+numtocombine],x_values_not_averaged[i],tolerance):
                currentList.append(x_values_not_averaged[i+numtocombine])
                numtocombine+=1
            i+=numtocombine
            x_values.append(self.avgfn([currentList]))
        if i == len(x_values_not_averaged)-1:
            x_values.append(x_values_not_averaged[-1]) 
        return np.array(x_values)

    def process_errbars_onlyvertical(self,xarys,yarys,tolerance):
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
                yvals = yarys[i][self.values_within_tolerance(xarys[i],x_value,tolerance)]
                if len(yvals) == 0:
                    mask[j] = False
                    continue
                ys[i][j] = self.avgfn(yvals)
                yerrs[i][j] = self.errfn(yvals)
            xs[i] = xs[i][mask]
            ys[i] = ys[i][mask]
            yerrs[i] = yerrs[i][mask]
        return xs,ys,yerrs

    def process_errbars_vertandhoriz(self,xarys,yarys):
        xs = np.empty(len(yarys),dtype = object)
        ys = np.empty(len(yarys),dtype = object)
        xerrs = np.empty(len(yarys),dtype = object)
        yerrs = np.empty(len(yarys),dtype = object)
        for i in range(len(xs)):
            xs[i] = np.zeros(len(xarys[i]))
            ys[i] = np.zeros(len(xarys[i]))
            xerrs[i] = np.zeros((len(xarys[i]),2))
            yerrs[i] = np.zeros((len(xarys[i]),2))
            for j,xlist in enumerate(xarys[i]):
                xs[i][j] = self.avgfn(xlist)
                xerrs[i][j] = self.errfn(xlist)
            for j,ylist in enumerate(yarys[i]):
                ys[i][j] = self.avgfn(ylist)
                yerrs[i][j] = self.errfn(ylist)
        return xs,ys,xerrs,yerrs

    def get_savefig_name(self,ion,labels,xVar,plot_type,custom_name = None):
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
        elif plot_type in [4]:
            # not implemented yet, just pretend it's type 1
            return self.get_savefig_name(ion,labels,xVar,1,custom_name)
        elif plot_type in [5]:
            ionNameNoSpaces = ion.replace(" ","")
            name = "%s_HeatMap_%s_%s" % (self.currentQuasarArrayName, ionNameNoSpaces, xVar)
        if name[0] == "_":
            name = name[1:]
        if custom_name:
            name = custom_name
        return name
    
    def get_redundancy_from_labels(self,labels):
        def allstartwith(ls,st):
            for s in ls:
                if not s.startswith(st):
                    return False
            return True
        
        numsplits = len(labels[0].split(":"))
        beginning = ""
        for i in range(numsplits):
            if not allstartwith(labels,beginning+labels[0].split(":")[i]):
                break
            else:
                beginning += labels[0].split(":")[i]
                beginning += ":"
        retlabels = [None]*len(labels)
        for i in range(len(retlabels)):
            retlabels[i] = labels[i].replace(beginning,"")
        if beginning.endswith(":"):
            beginning = beginning[:-1]
        return retlabels,beginning

    def get_ylabel_cd(self,ion,islogy):
        if ion in intensives:
            ylabel = intensiveslabels[ion]
            cd = ""
        elif ion in param_xVars:
            ylabel = param_unit_labels[yVar]
            cd = ""
        elif ":" in ion and ion.split(":")[1] != "cdens":
            ylabel = "%sfraction of %s in state"%(islogy,ion.split(":")[0])
            cd = "fraction "
        else:
            ylabel = "%scol dens"%islogy
            cd = "Column Density "
        return ylabel,cd  

    def plot_err(self, ion, quasarArray = None, xVar = "rdivR", save_fig = False, \
                 labels = None,extra_title = "",rlims = None,tolerance = 1e-5, \
                 dots = False,logx = False,logy = True, average = None,custom_name = None, \
                 coloration = None, visibility_threshold = None):
        print("Current constraints (name): "+self.currentQuasarArrayName)
        if not average is None:
            oldplots = self.plots
            self.setPlots(average)
        plt.figure()
        if isinstance(ion,list):
            ions = ion
            plot_type = 0
            for ion in ions:
                assert not ion in intensives
            assert quasarArray is None
            assert labels is None
            labels = ions
            xerr = None
            xarys,yarys = self.get_xy_type0(xVar,ions,quasarArray,rlims)
        elif not quasarArray is None and xVar in sightline_xVars:
            plot_type = 1
            assert isinstance(ion,str)
            assert len(labels) == len(quasarArray)
            xerr = None
            xarys,yarys = self.get_xy_type1(xVar,ion,quasarArray,rlims)
        elif (not quasarArray is None) and xVar in param_xVars and (not ion in param_xVars):
            plot_type = 2
            assert isinstance(ion,str)
            assert len(labels) == len(quasarArray)
            xerr = None
            xarys,yarys = self.get_xy_type2(xVar,ion,quasarArray,rlims)
        elif ion in param_xVars and xVar in param_xVars:
            yVar = ion
            plot_type = 3
            assert not quasarArray is None
            assert len(labels) == len(quasarArray)
            if not dots is None: 
                dots = True
            xerr = None
            xarys,yarys = self.get_xy_type3(xVar,yVar,quasarArray)
        elif xVar not in sightline_xVars + param_xVars:
            yVar = ion
            plot_type = 4
            xarys,yarys = self.get_xy_type4(xVar,yVar,quasarArray,rlims)
        else:
            print("Requirements not met (type could not be detected). Cannot plot.")
            return
        
        labels, redundancy = self.get_redundancy_from_labels(labels)

        if plot_type in [0,1,2,3]:
            for i in range(len(yarys)):
                xarys[i] = xarys[i][yarys[i]>0]
                yarys[i] = yarys[i][yarys[i]>0]
                if logy:
                    yarys[i] = np.log10(yarys[i])
        elif plot_type in [4]:
            for i in range(len(yarys)):
                for j in range(len(yarys[i])):
                    xarys[i][j] = xarys[i][j][yarys[i][j]>0]
                    yarys[i][j] = yarys[i][j][yarys[i][j]>0]
                    yarys[i][j] = yarys[i][j][xarys[i][j]>0]
                    xarys[i][j] = xarys[i][j][xarys[i][j]>0]
                    if logy:
                        yarys[i][j] = np.log10(yarys[i][j])
                    if logx:
                        xarys[i][j] = np.log10(xarys[i][j])
        if visibility_threshold:
            if plot_type in [0,1,2,3]:
                xarys_visible = np.empty(len(xarys),dtype = object)
                yarys_visible = np.empty(len(xarys),dtype = object)
                for i in range(len(yarys)):
                    oldyarys = np.power(10,yarys[i]) if logy else yarys[i]
                    xarys_visible[i] = xarys[i][oldyarys>visibility_threshold]
                    yarys_visible[i] = yarys[i][oldyarys>visibility_threshold]
            elif plot_type in [4]:
                xarys_visible = np.empty(len(xarys),dtype = object)
                yarys_visible = np.empty(len(xarys),dtype = object)
                for i in range(len(yarys)):
                    xarys_visible[i] = np.empty(len(yarys[i]), dtype = object)
                    yarys_visible[i] = np.empty(len(yarys[i]), dtype = object)
                    for j in range(len(yarys[i])):
                        oldyarys = np.power(10,yarys[i][j]) if logy else yarys[i][j]
                        oldxarys = np.power(10,xarys[i][j]) if logx else xarys[i][j]
                        xarys_visible[i][j] = xarys[i][j][oldyarys>visibility_threshold[0]]
                        yarys_visible[i][j] = yarys[i][j][oldyarys>visibility_threshold[0]]
                        yarys_visible[i][j] = yarys[i][j][oldxarys>visibility_threshold[1]]
                        xarys_visible[i][j] = xarys[i][j][oldxarys>visibility_threshold[1]]
        if not average is None:
            oldplots = self.plots
            self.setPlots(average)
        if plot_type in [0,1,2,3]:
            xs,ys,yerrs = self.process_errbars_onlyvertical(xarys,yarys,tolerance)
            xerrs = None
            if visibility_threshold:
                xs_visible,ys_visible,yerrs_visible = self.process_errbars_onlyvertical(xarys_visible,yarys_visible,tolerance)
                xerrs_visible = None
        elif plot_type in [4]:
            xs,ys,xerrs,yerrs = self.process_errbars_vertandhoriz(xarys,yarys)
            if visibility_threshold:
                xs_visible,ys_visible,xerrs_visible,yerrs_visible = self.process_errbars_vertandhoriz(xarys_visible,yarys_visible)
        islogy = "log " if logy else ""
        ylabel, cd = self.get_ylabel_cd(ion,islogy)
        if plot_type in [1,2,3,4]:
            ylabel = ion +" "+ ylabel
        plt.ylabel(ylabel)
        if xVar in sightline_xVars:
            plt.xlabel(sightline_unit_labels[xVar])
        elif xVar in param_xVars:
            plt.xlabel(param_unit_labels[xVar])
        if plot_type in [0]:
            ionstr = redundancy
        elif plot_type in [1,2,3]:
            ionstr = ion
        elif plot_type in [4]:
            ionstr = ""
        ionstr = ion if plot_type in [1,2,3] else redundancy
        plt.title('%s %sAverages (%s) %s'%(ionstr, cd, self.plots, extra_title))


        for i in range(len(xs)):
            x = xs[i]
            y = ys[i]
            if plot_type == 4:
                islogx = "log " if logx else ""
                if xVar in sightline_xVars:
                    plt.xlabel(sightline_unit_labels[xVar])
                else:
                    xlabel, _ = self.get_ylabel_cd(xVar,islogx)
                    plt.xlabel(xVar +" "+ xlabel)
            yerr = np.transpose(yerrs[i])
            xerr = np.transpose(xerrs[i]) if not (xerrs is None) else None
            label = labels[i] if not visibility_threshold else None
            color = None
            if not (coloration is None):
                color = coloration[i]
            #change fmt to . or _ for a dot or a horizontal line
            fmtdict = {"mean":',',"median_std":',',"median":"."}
            alpha_value = .3 if visibility_threshold else 1.
            if self.plots == "scatter" and dots:
                print "average = scatter AND dots = True detected"
                print "'scatter' gives all sightlines without averaging."
                print "'dots' gives all averages without error bars."
                print "please change one and try again"
                return
            elif dots:
                plt.plot(x,y,"o",label = label,color = color,alpha = alpha_value)
                if visibility_threshold:
                    x_visible = xs_visible[i]
                    y_visible = ys_visible[i]
                    plt.plot(x_visible,y_visible,"o",color = color,label = labels[i])
            elif self.plots == "scatter":
                def flatten_if_needed(ary):
                    try:
                        for item in ary:
                            len(item)
                        return np.concatenate(ary)
                    except:
                        return ary
                x = flatten_if_needed(xarys[i])
                y = flatten_if_needed(yarys[i])
                plt.plot(x,y,'o',label = label, markersize = 2,color = color,alpha = alpha_value)
                if visibility_threshold:
                    x_visible = flatten_if_needed(xarys_visible[i])
                    y_visible = flatten_if_needed(yarys_visible[i])
                    plt.plot(x_visible,y_visible,"o",label = labels[i],markersize = 2,color = color)
            else:
                plt.errorbar(x,y, xerr = xerr, yerr=yerr, fmt=fmtdict[self.plots], capsize = 3, label = label,color = color,alpha = alpha_value)
                if visibility_threshold:
                    x_visible = xs_visible[i]
                    y_visible = ys_visible[i]
                    yerr_visible = np.transpose(yerrs_visible[i])
                    xerr_visible = np.transpose(xerrs_visible[i]) if not (xerrs_visible is None) else None
                    plt.errorbar(x_visible,y_visible, xerr = xerr_visible, yerr=yerr_visible, fmt=fmtdict[self.plots], capsize = 3, label = labels[i],color = color)
        plt.legend()
        
        if logx and not plot_type == 4:
            plt.xscale('log')
            
        if save_fig:
            name = self.get_savefig_name(ion,labels,xVar,plot_type,custom_name = custom_name)
            plt.savefig("plots/"+name + ".png")
            return plt,"plots/"+name + ".png"
        
        if not average is None:
            self.setPlots(oldplots)
        
        return plt

    def definecolorbar(self):
        from matplotlib.colors import LinearSegmentedColormap
        f= 256.0
        cdict = {'red':   ((0.0,  255/f, 255/f),
                           (1e-9, 255/f, 255/f),
                           (0.1,  255/f, 255/f),
                           (0.3,  255/f, 255/f),
                           (0.5,  139/f, 139/f),
                           (1.0,  0/f, 0/f)),

                 'green': ((0.0,  255/f, 255/f),
                           (1e-9, 255/f, 255/f),
                           (0.1,  165/f, 165/f),
                           (0.3,  0/f, 0/f),
                           (0.5,  0/f, 0/f),
                           (1.0,  0/f, 0/f)),

                 'blue':  ((0.0,  255/f, 255/f),
                           (1e-9, 255/f, 0/f),
                           (0.1,  0/f, 0/f),
                           (0.3,  0/f, 0/f),
                           (0.5,  0/f, 0/f),
                           (1.0,  0/f, 0/f))}

        hotcustom = LinearSegmentedColormap('HotCustom', cdict)
        return hotcustom

    #summary: plots histogram(s) of ion column density vs incrementing xvariable, 
    #         uses a color bar to show the percentage of sightlines for a certain column density at a specific x value
    #
    #    def plot_err(self, ion, quasarArray = None, xVar = "r", save_fig = False, \
    #             reset = False, labels = None,extra_title = "",rlims = None,\
    #             tolerance = 1e-5,dots = False,logx = False,average = None,logy = True
    def plot_hist(self, ion, xVar = "rdivR",extra_title = "",rlims = None,weights = True,\
                  save_fig = False, tolerance = 1e-5,ns = (42,15),logx = False, logy = True):
        hotcustom = self.definecolorbar()
        plt.register_cmap(cmap=hotcustom)
        if rlims is None:
            rlims = [0.1,np.inf]
        vardict = {"theta":1,"phi":2,"r":3,"rdivR":3}
        distances = "kpc" if xVar == "r" else "Rvir"
        gq = GeneralizedQuasarSphere(self.currentQuasarArray,distance = distances)
        xs = gq.info[:, vardict[xVar]]
        ys = self.get_yVar_from_str(gq,ion)
        rs = gq.info[:,3]
        acceptedLines = np.logical_and(rlims[0]<=rs,rs<=rlims[1])
        xs = xs[acceptedLines]
        ys = ys[acceptedLines]
        islogy = ""
        if logy:
            xs = xs[ys>0]
            ys = ys[ys>0]
            ys = np.log10(ys)
            islogy = "log "
        if weights:
            weight = xs*0.0
            for i in range(len(xs)):
                if xs[i]!=0:
                    weight[i] = 1.0/len(xs[abs(xs[i]/xs - 1)<tolerance])
                else:
                    weight[i] = 1.0/len(xs[xs<tolerance])
            cbarlabel = "Fraction of lines for fixed %s"%(xVar)
        else:
            weight = xs*0.0+1.0
            cbarlabel = "Total number of lines"
        H, xedges, yedges = np.histogram2d(xs, ys, bins=ns,weights = weight)
        H = H.T
        X, Y = np.meshgrid(xedges, yedges)
        plt.pcolormesh(X,Y, H, cmap=hotcustom)
        plt.colorbar(label = cbarlabel)
        ylabel,cd = self.get_ylabel_cd(ion,islogy)
        plt.ylabel(ylabel)
        plt.xlabel(sightline_unit_labels[xVar])
        plt.title('%s %sDistribution %s'%(ion, cd, extra_title))
        def get_new_axislim(current,divideBy = 20):
            x1 = current[0]
            x2 = current[1]
            return x1-(x2-x1)/divideBy,x2+(x2-x1)/divideBy
        plt.ylim(get_new_axislim(plt.ylim()))
        if logx:
            plt.xscale('log')
            plt.xlabel('log ' + sightline_unit_labels[xVar])
        else:
            plt.xlim(get_new_axislim(plt.xlim()))

        return plt

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

def print_general_shape(ary,name = "ary",indent=0,showitems = False):
    if indent == 0:
        print name+": len()="+str(len(ary))
        indent += 1
    for i,a in enumerate(ary):
        try:
            print '---_'*indent+str(len(a))
            print_general_shape(a,indent = indent+1)
        except:
            if showitems:
                print '---_'*indent+"*item*"+str(a)
