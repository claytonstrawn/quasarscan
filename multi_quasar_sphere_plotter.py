import numpy as np
import os
import sys
import matplotlib.pyplot as plt
import datetime
from functools import reduce
try:
    from quasarscan import quasar_sphere
    from quasarscan import observational_quasar_sphere
    from quasarscan import ion_lists
    from quasarscan import gasbinning
    from quasarscan import roman
    from quasarscan.sorter import MultiSphereSorter
    level = 0
except:
    import quasar_sphere
    import observational_quasar_sphere
    import ion_lists
    import gasbinning
    import roman
    from sorter import MultiSphereSorter
    level = 1
#precondition: assumes there are only two levels of depth within the output folder
#postcondition: returns a list of all textfiles
def get_all_textfiles(inquasarscan = 1,loadonly = 'all'):
    
    #pathname by default starts in output
    path = "output"
    if not inquasarscan:
        path = "quasarscan/output"
    textfiles = []
    def one_is_in_name(name,loadonly):
        if loadonly == 'all':
            return True
        if isinstance(loadonly,str):
            loadonly = [loadonly]
        for sim in loadonly:
            if sim in name:
                return True
        return False
    #gets all folders in output
    dirs = os.listdir(path)
    for folderName in dirs:
        if not folderName.startswith(".") and one_is_in_name(folderName,loadonly):
            folderPath = path + "/" + folderName
            folderDirs = os.listdir(folderPath)
            for fileName in folderDirs:
                if not fileName.startswith("."):
                    textfiles.append(os.path.join(folderPath,fileName))
    return textfiles

def get_all_observations(inquasarscan = 1,loadonly = 'all'):
    path = "observations"
    if not inquasarscan:
        path = "quasarscan/observations"
    textfiles = []
    def one_is_in_name(name,loadonly):
        if loadonly == 'all':
            return True
        if isinstance(loadonly,str):
            loadonly = [loadonly]
        for sim in loadonly:
            if sim in name:
                return True
        return False
    #gets all folders in output
    dirs = os.listdir(path)
    for fileName in dirs:
        if not fileName.startswith(".") and one_is_in_name(fileName,loadonly):
            textfiles.append(os.path.join(path,fileName))
    return textfiles
    
def sort_ions(ions,flat = True):
    def sort_ions_one_element(ions,element):
        nums = [None]*len(ions)
        toreturn = []
        for i in range(len(ions)):
            nums[i] = roman.from_roman(ions[i].split(" ")[1])
        nums.sort()
        for val in nums:
            toreturn.append("%s %s"%(element,roman.to_roman(val)))
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

def reversearray(ary):
    ary = list(ary)
    ary.reverse()
    return np.array(ary)

stringcriteria = ["ions","name","simname","version","code","simnum","Rvir_is_real","compaction_stage"]
intensives = ["Z","T","rho"]
intensiveslabels = {"Z":"avg metallicity","T":"avg temperature","rho":"avg density"}
intensivespositions = {"Z":-1,"T":-2,"rho":-3}
sightline_xVars = ["r","rdivR","theta","phi"]
param_xVars = ["redshift","a0","Mvir","gas_Rvir","star_Rvir","dm_Rvir","sfr","ssfr","L_mag","Mstar","Mgas","Rvir"]
sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
           "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)","phi" \
           :"azimuthal viewing angle (rad)"}
param_unit_labels = {"redshift":"z","a0":"a","Rvir":'Virial radius (kpc)',"Mvir":"Virial Mass (Msun)",\
                    "gas_Rvir":"Gas Mass within Rvir (Msun)","Mgas":"Gas Mass within Rvir (Msun)","star_Rvir":"Stellar Mass within Rvir (Msun)",\
                    "Mstar":"Stellar Mass within Rvir (Msun)","dm_Rvir":"Dark Matter Mass within Rvir (Msun)","sfr":"Star Formation Rate (Msun yr-1)",\
                    "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)","L_mag":"Magnitude of Angular Momentum"}

class MultiQuasarSpherePlotter():
    #param: textfiles     if a list of textfiles is specified, those specific textfiles will be loaded; else,
    #                     all textfiles in output are loaded
    def __init__(self, loadonly = "all",loadobs = 'all',textfiles = None, cleanup = False,plots = "mean",throwErrors = False,safetycheck = True):
        self.plots = "mean"
        self.avgfn = np.mean
        self.setPlots(plots)
        self.quasarArray = []
        if textfiles is None:
            textfiles = get_all_textfiles(level,loadonly = loadonly)
            observations = get_all_observations(level,loadonly = loadobs)
        for textfile in textfiles:
            try:
                readvalsoutput = quasar_sphere.read_values(textfile)
                q = quasar_sphere.QuasarSphere(readvalsoutput = readvalsoutput)
                if safetycheck and self.pass_safety_check(q):
                    self.quasarArray.append(q)
                elif cleanup:
                    todo = input("file %s did not pass safety check. Remove it? (y/n)"%textfile).lower()
                    os.remove(textfile) if todo == 'y' else None
            except Exception as e:
                print(textfile + " could not load because:")
                print(e)
                if throwErrors:
                    raise e
        self.observationArray = np.array([])
        for textfile in observations:
            papername,headers,lines = observational_quasar_sphere.read_textfile(textfile)
            self.observationArray = np.append(self.observationArray,\
                      observational_quasar_sphere.create_observational_quasarArray(papername,headers,lines))
                                
        self.quasarArray = np.array(self.quasarArray)
        self.currentObservationArray = np.copy(self.observationArray)
        self.currentQuasarArray = np.copy(self.quasarArray)
        self.currentQuasarArrayName = ''
        
        if len(self.currentQuasarArray) == 0:
            print ("There are no quasarspheres stored in currentQuasarArray!")
    
    def length(self):
        return len(self.currentQuasarArray)
    def list_all_QuasarSpheres(self, *criteria):
        s = ""
        if len(self.currentQuasarArray)==0:
            print("no QuasarSpheres")
        for q in sorted(self.currentQuasarArray, key=lambda x: (x.name,x.rounded_redshift)):
            s += q.name +" at z=%s ("%q.rounded_redshift
            for c in criteria:
                v = eval('q.%s'%c)
                if isinstance(v,str) or isinstance(v,list):
                    s+="%s = %s, "%(c,v)
                elif .001<v<1:
                    s+="%s = %0.3f, "%(c,v)
                elif 1<v<1000:
                    s+="%s = %3.0f, "%(c,v)
                else:
                    s+="%s = %1.3e, "%(c,v)
            s = s[:-2]
            if len(criteria)>0:
                s+=")"
            s+='\n'
        print(s[:-1])
            
    #tests each condition, each condition is a safety check
    def pass_safety_check(self, q,require_metadata = False):
        minlength = 10
        minions = []
        if q.length_reached < minlength:
            print("Length for %s at redshift %s is not valid." %(q.name, str(q.rounded_redshift)))
            return False
        elif not all(x in q.ions for x in minions):
            print("Not all necessary ions present in %s at redshift %s."%(q.name, str(q.rounded_redshift)))
            return False
        if require_metadata:
            if q.final_a0 is None:
                print("metadata for %s at redshift %s is not valid." %(q.name,str(q.rounded_redshift)))
                return False
        return True
        
    def reset_current_Quasar_Array(self):
        self.currentQuasarArray = np.copy(self.quasarArray)
        self.currentObservationArray = np.copy(self.observationArray)
        self.currentQuasarArrayName = ''
        
    def check_criteria_works(self,constrainCriteria,atEnd = False,**kwargs):
        if len(self.currentQuasarArray)==0:
            #"Cannot constrain further"
            return
        elif constrainCriteria not in self.currentQuasarArray[0].__dict__.keys():
            print("Constrain criteria " + constrainCriteria + " does not exist. Please re-enter a valid criteria.")
            raise Exception
        elif isinstance(atEnd,float):
            oldQuasarArray=np.copy(self.currentQuasarArray)
            self.constrain_current_Quasar_Array("final_a0",bins=[atEnd-.1,np.inf],changeArrayName=False)
            if len(self.currentQuasarArray) == 0:
                print("No galaxies get to that high of a0")
            return oldQuasarArray
        
    def prepare_to_sort(self, criteria, bins,atEnd,**kwargs):
        
        oldQuasarArray = self.check_criteria_works(criteria,atEnd=atEnd,**kwargs)
        if criteria == "ions":
            raise Exception("You cannot sort by 'ions'")
            
        # can sort by giving one value, it will assume you mean "less than value" and "greater than value"
        if isinstance(bins, float) or isinstance(bins, int):
            bins = np.array([0.0, bins, np.inf])
        elif isinstance(bins, str) and criteria in stringcriteria:
            bins = [bins]
        
        return bins,oldQuasarArray
    
    def postprocess_sorted(self, labels, bins, quasarBins, onlyNonempty = False,reverse = False,**kwargs):
        if quasarBins is None:
            raise Exception("No quasars in quasarBin!")
        empty = True
        nonemptyArray = []
        nonemptyLabelArray = []
        for i,item in enumerate(quasarBins):
            if len(item)>0:
                nonemptyArray.append(item)
                nonemptyLabelArray.append(labels[i])
                empty = False
        if onlyNonempty:
            labels, bins, quasarBins = np.array(nonemptyLabelArray),bins, np.array(nonemptyArray)
            print("All bins are empty." if empty else "")
        if reverse:
            labels = reversearray(labels)
            bins = reversearray(bins)
            quasarBins = reversearray(quasarBins)
            
        return labels, bins, quasarBins
    
    
    #param bins either serves as an array with each element being as a cutpoint, or a single value
    #follows array indexing exclusivity and inclusivity rules
    def sort_by(self, criteria, bins = [0,np.inf],\
        atEnd = False,splitEven = 0,**kwargs):
        bins,oldQuasarArray = self.prepare_to_sort(criteria,bins,atEnd,**kwargs)
        sorter = MultiSphereSorter(self.currentQuasarArray)
        obs_sorter = MultiSphereSorter(self.currentObservationArray)
        if splitEven:
            labels, bins, quasarBins = sorter.splitEven(criteria,splitEven,atEnd = atEnd)
        else:
            labels, bins, quasarBins = sorter.sort(criteria,bins,atEnd = atEnd)

        _, _, obsBins = obs_sorter.sort(criteria,bins,atEnd = False)
        labels,bins,quasarBins = self.postprocess_sorted(labels,bins,quasarBins,**kwargs)
        _,_,obsBins = self.postprocess_sorted(labels,bins,obsBins,**kwargs)
        return labels,bins, quasarBins, obsBins
    
    def constrain_via_gasbins(self,gasbintype=None):
        if gasbintype == None:
            gasbintype = input("Available bins are: %s"%gasbinning.possible_bin_types)
        g = gasbinning.GasBinsHolder(bins=[gasbintype])
        toReturn = []
        for q in self.currentQuasarArray:
            if g.get_bin_str() in q.gasbins.get_bin_str():
                toReturn.append(q)
        self.currentQuasarArray = np.array(toReturn)
        
    def get_bin_values(self,constrainCriteria,bins,exclude=False,**kwargs):
        if isinstance(bins, list):                
            if len(bins) != 2 and constrainCriteria not in stringcriteria:
                print("Length of bins must be 2: [lower,upper]")
                raise Exception
        elif isinstance(bins,str) and constrainCriteria in stringcriteria:
            bins = [bins]
        if exclude:
            possible_bins = []
            for q in self.currentQuasarArray:
                possible_bins.append(eval("q.%s"%constrainCriteria))
            possible_bins = set(possible_bins)
            excluded_bins = bins
            bins = list(possible_bins.difference(set(bins)))
        return bins
    
    def constrain_array_helper(self,sorter,obs_sorter,constrainCriteria,bins,splitEven=False,atEnd=False,\
                               set_main_array=False,sortobs='default',**kwargs):
        if not splitEven:
            labels, bins, temp = sorter.sort(constrainCriteria,bins,atEnd=atEnd)
            _,_,obs_temp = obs_sorter.sort(constrainCriteria,bins)
        else:
            labels,bins,temp = sorter.splitEven(constrainCriteria,2,atEnd=atEnd)
            _,_,obs_temp = obs_sorter.sort(constrainCriteria,bins)
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
            obs_temp = np.array([obs_temp[take]])
        self.currentQuasarArray = np.concatenate(temp)
        if sortobs == True or (sortobs == 'default' and (constrainCriteria not in stringcriteria or constrainCriteria == 'ions')):
            self.currentObservationArray = np.concatenate(obs_temp)
        if set_main_array:
            self.quasarArray = np.copy(self.currentQuasarArray)
        return bins
    
    def change_array_name(self,constrainCriteria,bins,changeArrayName=True,exclude=False,**kwargs):
         if changeArrayName:
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
                if lowlabel == "lessthan" and highlabel == "andhigher":
                    self.currentQuasarArrayName += "any"
                else:
                    self.currentQuasarArrayName += "%s%s"%(lowlabel,highlabel)
            else:
                if exclude:
                    self.currentQuasarArrayName += 'exclude'
                for acceptedValue in bins:
                    self.currentQuasarArrayName += acceptedValue.replace(" ","")

    def constrain_current_Quasar_Array(self, constrainCriteria,bins=None,**kwargs):
        self.check_criteria_works(constrainCriteria,**kwargs)
        bins = self.get_bin_values(constrainCriteria,bins,**kwargs)
        sorter = MultiSphereSorter(self.currentQuasarArray)
        obs_sorter = MultiSphereSorter(self.currentObservationArray)
        bins = self.constrain_array_helper(sorter,obs_sorter,constrainCriteria,bins,**kwargs)
        self.change_array_name(constrainCriteria,bins,**kwargs)
        return bins

    def setPlots(self,plots,quartiles = None,**kwargs):
        if isinstance(plots,tuple):
            if len(plots) == 3:
                quartiles = (plots[1],plots[2])
            elif len(plots) == 2:
                if isinstance(plots[1],tuple):
                    quartiles = plots[1]
                elif plots[0] == 'median' or plots[0] == 'med':
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
        def getstddev(data):
            return np.array([np.nanstd(data),np.nanstd(data)])
        if plots in ["mean","stderr"]:
            self.plots = "mean"
            self.avgfn = np.nanmean
            self.errfn = getstderr
        if plots in ["stddev"]:
            self.plots = "stddev"
            self.avgfn = np.nanmean
            self.errfn = getstddev
        elif plots in ["median","med"]:
            self.plots = "median"
            self.avgfn = np.nanmedian
            self.errfn = getquartiles
        elif plots in ["med_noquartiles","median_std"]:
            self.plots = "median_std"
            self.avgfn = np.nanmedian
            self.errfn = getstderr
        elif plots in ["scatter"]:
            self.plots = "scatter"
            self.avgfn = np.nanmean
            self.errfn = getstderr
        elif plots in ['covering_fraction','cvf']:
            self.plots = "covering_fraction"
            self.avgfn = np.nanmean
            self.errfn = getstderr
            
    def get_yVar_from_str(self,gq,stringVar):
        def split_by_ops(s):
            #NOTE: This will have a bug if a variable legitimately ends in e and adds/subtracts
            #the next term! this is because e-XX is how python floats are written as strings
            s = s.replace("e-","_temp_minus_char_")
            s = s.replace("e+","_temp_plus_char_")
            s = s.replace("(","_splitchar_")
            s = s.replace(")","_splitchar_")
            s = s.replace("/","_splitchar_")
            s = s.replace("*","_splitchar_")
            s = s.replace("+","_splitchar_")
            s = s.replace("-","_splitchar_")
            s = s.replace("_temp_minus_char_","e-")
            s = s.replace("_temp_plus_char_","e+")
            return filter(None,s.split("_splitchar_"))
        strings_to_find = split_by_ops(stringVar)
        new_str_to_eval = stringVar
        strings_to_replace_with = {}
        for i,s in enumerate(strings_to_find):
            if len(s)>0:
                try:
                    eval(s)
                    string_to_replace_with = s
                except:
                    try:
                        string_to_replace_with = "gq.info[:,%d]"%gq.get_ion_column_num(s)
                    except:
                        continue
                strings_to_replace_with[s] = string_to_replace_with
        for s in sorted(strings_to_replace_with.keys(),key=len,reverse=True):
            new_str_to_eval = new_str_to_eval.replace(s,strings_to_replace_with[s])
        to_return = eval(new_str_to_eval)
        return to_return

    def get_xy_type0(self,xVar,yVars,rlims):
        quasarArray = [self.currentQuasarArray]
        xs = np.empty(len(yVars),dtype = object)
        ys = np.empty(len(yVars),dtype = object)
        for i in range(len(yVars)):
            yVar = yVars[i]
            if xVar in sightline_xVars:
                x,y = self.get_xy_type1(xVar,yVar,quasarArray,rlims)
            elif xVar in param_xVars:
                x,y = self.get_xy_type2(xVar,yVar,quasarArray,rlims)
            else:
                print("Not type 2!")
            xs[i] = x[0]
            ys[i] = y[0]
        return xs,ys

    def get_xy_type1(self,xVar,yVar,quasarArray,rlims):
        if rlims is None:
            rlims = [0.1,np.inf]
        elif rlims == "all":
            rlims = [0.0,np.inf]
        vardict = {"theta":1,"phi":2,"r":3,"rdivR":3}
        distances = "kpc" if xVar == "r" else "Rvir"
        gqary = []
        xs = np.empty(len(quasarArray),dtype = object)
        ys = np.empty(len(quasarArray),dtype = object)
        for i in range(len(quasarArray)):
            gq = quasar_sphere.GeneralizedQuasarSphere(quasarArray[i],distance=distances)
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
        elif rlims == "all":
            rlims = [0.0,np.inf]
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
                acceptedLines = np.logical_and(rlims[0]<=rs/q.Rvir,\
                                               rs/q.Rvir<=rlims[1])
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
        elif rlims is "all":
            rlims = [0.0,np.inf]
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
                acceptedLines = np.logical_and(rlims[0]<=rs/q.Rvir,\
                                               rs/q.Rvir<=rlims[1])
                x = x[acceptedLines]
                y = y[acceptedLines]
                xs[i][j] = x
                ys[i][j] = y
        return xs,ys

    def values_within_tolerance(self,x_in_list,x_comp,tolerance,symmetric = True):
        tocompare = x_in_list-x_comp
        if symmetric:
            tocompare = np.abs(tocompare)
        return np.logical_and(-tolerance<=tocompare,tocompare<=tolerance)

    def combine_xs(self,x_variable,tolerance):
        if tolerance == 0:
            return x_variable
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
            x_values.append(np.min([currentList]))
        if i == len(x_values_not_averaged)-1:
            x_values.append(x_values_not_averaged[-1]) 
        return np.array(x_values)

    def decide_plot_type(self,ion,xVar,labels=None,quasarArray=None,lq=None,**kwargs):
        if lq:
            labels = lq[0]
        if isinstance(ion,tuple):
            ion = ion[0]
        if isinstance(ion,list):
            assert labels is None and quasarArray is None
            return 0
        elif isinstance(ion,str) and xVar in sightline_xVars:
            return 1
        elif xVar in param_xVars and (not ion in param_xVars):
            assert isinstance(ion,str)
            return 2
        elif ion in param_xVars and xVar in param_xVars:
            return 3
        elif xVar not in sightline_xVars + param_xVars:
            return 4
        else:
            print("Requirements not met (type could not be detected). Cannot plot.")

    def get_labels_from_ion(self,plot_type,ion,xVar,labels=None,lq=None,**kwargs):
        if lq is not None:
            labels=lq[0]
        if isinstance(ion,tuple):
            ion_name=ion[1]
            ion = ion[0]
        else:
            ion_name=ion
        if isinstance(xVar,tuple):
            xVar_name = xVar[1]
            xVar = xVar[0]
        else:
            xVar_name = xVar
        if plot_type==0:
            ions = ion
            labels = []
            ions_notuple = []
            for ion in ions:
                if isinstance(ion,tuple):
                    ion_name = ion[1]
                    ion = ion[0]
                else:
                    ion_name = ion
                labels.append(ion_name)
                ions_notuple.append(ion)
                assert not ion in intensives
            ion = ions_notuple
            ion_name = labels
        if labels is None:
            labels = ['all']
        return ion,ion_name,xVar,xVar_name,labels

    def get_sightline_xy_vals(self,plot_type,ion,xVar,lq=None,quasarArray=None,rlims=None,**kwargs):
        if lq is not None:
            quasarArray = lq[2]
        if quasarArray is None:
            quasarArray = [self.currentQuasarArray]
        if rlims is None:
            rlims = np.array([0.1,np.inf])
        elif rlims == "all":
            rlims = [0.0,np.inf]
        if xVar == 'rdivR':
            self.constrain_current_Quasar_Array("Rvir_is_real",bins=['True'],changeArrayName=False)
        if plot_type==0:
            xarys,yarys = self.get_xy_type0(xVar,ion,rlims)
        elif plot_type==1:
            xarys,yarys = self.get_xy_type1(xVar,ion,quasarArray,rlims)
        elif plot_type==2:
            xarys,yarys = self.get_xy_type2(xVar,ion,quasarArray,rlims)
        elif plot_type==3:
            xarys,yarys = self.get_xy_type3(xVar,ion,quasarArray)
        elif plot_type==4:
            xarys,yarys = self.get_xy_type4(xVar,ion,quasarArray,rlims)

        if plot_type in [0,1,2,3]:
            for i in range(len(yarys)):
                xarys[i] = xarys[i][np.logical_and(yarys[i]>=0,yarys[i]<np.inf)]
                yarys[i] = yarys[i][np.logical_and(yarys[i]>=0,yarys[i]<np.inf)]
        elif plot_type in [4]:
            for i in range(len(yarys)):
                for j in range(len(yarys[i])):
                    xarys[i][j] = xarys[i][j][np.logical_and(yarys[i][j]>=0,yarys[i][j]<np.inf)]
                    yarys[i][j] = yarys[i][j][np.logical_and(yarys[i][j]>=0,yarys[i][j]<np.inf)]
                    yarys[i][j] = yarys[i][j][np.logical_and(xarys[i][j]>=0,xarys[i][j]<np.inf)]
                    xarys[i][j] = xarys[i][j][np.logical_and(xarys[i][j]>=0,xarys[i][j]<np.inf)]
                    if xVar == 'rho':
                        xarys[i][j]/=1.67373522e-24
        self.debug = xarys,yarys
            
        return xarys,yarys
    
    def should_take_logs_xy(self,ion,xVar,logx,logy,average='default',**kwargs):
        probablylinear = ['redshift','rounded_redshift','a0','Rvir']
        if logx=='guess':
            if xVar in sightline_xVars:
                logx=False
            elif xVar in probablylinear:
                logx=False
            elif xVar in param_xVars:
                logx=True
            else:
                logx=True
        if logy=='guess':
            if ion in sightline_xVars or average=='covering_fraction':
                logy=False
            elif ion in param_xVars and xVar not in probablylinear:
                logy=True
            else:
                logy=True
        return logx,logy
    

    def process_scatter_points(self,xarys,yarys,logx,logy,subsample=1.0,offsetx=False,tolerance=1e-5,**kwargs):
        def flatten_if_needed(ary):
            try:
                for item in ary:
                    len(item)
                return np.concatenate(ary)
            except:
                return ary
        def getsubsample(data,frac):
            l = int(len(data)*frac)
            mask = np.zeros(len(data))
            mask[:l] = 1
            np.random.shuffle(mask)
            return mask.astype(bool)
        xs = np.zeros(len(xarys),dtype=object)
        ys = np.zeros(len(xarys),dtype=object)
        empty = True
        for i in range(len(xarys)):
            xs[i] = flatten_if_needed(xarys[i])
            ys[i] = flatten_if_needed(yarys[i])
            mask = getsubsample(xs[i],subsample)
            xs[i]=xs[i][mask]
            ys[i]=ys[i][mask]
            if np.any(ys[i]>0) and np.any(xs[i]>0):
                empty=False
            if logx:
                xs[i] = np.log10(xs[i])
            if logy:
                ys[i] = np.log10(ys[i])
            if offsetx and np.any(ys[i]>0) and np.any(xs[i]>0):
                decimals = int(-np.log10(tolerance))
                randomscatterwidth = np.min(np.diff(np.unique(np.round(xs[i],decimals = decimals))))
                add = (np.random.random(len(xs[i]))-.5)*randomscatterwidth*.5
                xs[i] += add
        return xs,ys,empty

    def process_errbars_onlyvertical(self,xarys,yarys,logx,logy,offsetx=False,tolerance=1e-5,**kwargs):
        xs = np.empty(len(yarys),dtype = object)
        ys = np.empty(len(yarys),dtype = object)
        yerrs = np.empty(len(yarys),dtype = object)
        empty = True
        for i in range(len(yarys)):
            unique_xs = self.combine_xs(xarys[i],tolerance)
            xs[i] = unique_xs
            ys[i] = np.zeros(len(unique_xs))
            yerrs[i] = np.zeros((len(unique_xs),2))
            mask = np.ones(len(unique_xs),dtype = bool)
            for j in range(len(unique_xs)):
                x_value = unique_xs[j]
                yvals = yarys[i][self.values_within_tolerance(xarys[i],x_value,tolerance,symmetric=False)]
                if len(yvals) == 0:
                    mask[j] = False
                    continue
                ys[i][j] = self.avgfn(yvals)
                yerrs[i][j] = self.errfn(yvals)
                if ys[i][j] > 0:
                    empty = False
                if logy:
                    yerrs[i][j][0] = np.log10(ys[i][j])-np.log10(ys[i][j]-yerrs[i][j][0])
                    if yerrs[i][j][0] == np.inf or np.isnan(yerrs[i][j][0]):
                        yerrs[i][j][0] = max(np.log10(ys[i][j]),0)
                    yerrs[i][j][1] = np.log10(ys[i][j]+yerrs[i][j][1])-np.log10(ys[i][j])
                    ys[i][j] = np.log10(ys[i][j])
            xs[i] = xs[i][mask]
            ys[i] = ys[i][mask]
            yerrs[i] = yerrs[i][mask]
            if offsetx and len(xs[i])>0:
                scatterwidth = np.min(np.diff(np.unique(xs[i])))
                if logx:
                    xs[i][xs[i]>0]*=10**((float(i)/len(xs)-.5)*scatterwidth*.2)
                else:
                    xs[i]+=(float(i)/len(xs)-.5)*scatterwidth*.5
            if logx:
                xs[i] = np.log10(xs[i])
        return xs,ys,yerrs,empty

    def process_errbars_vertandhoriz(self,xarys,yarys,logx,logy,**kwargs):
        xs = np.empty(len(yarys),dtype = object)
        ys = np.empty(len(yarys),dtype = object)
        xerrs = np.empty(len(yarys),dtype = object)
        yerrs = np.empty(len(yarys),dtype = object)
        empty = True
        for i in range(len(xs)):
            xs[i] = np.zeros(len(xarys[i]))
            ys[i] = np.zeros(len(xarys[i]))
            xerrs[i] = np.zeros((len(xarys[i]),2))
            yerrs[i] = np.zeros((len(xarys[i]),2))
            for j,xlist in enumerate(xarys[i]):
                xs[i][j] = self.avgfn(xlist)
                xerrs[i][j] = self.errfn(xlist)
                if logx:
                    xerrs[i][j][0] = np.log10(xs[i][j])-np.log10(xs[i][j]-xerrs[i][j][0])
                    xerrs[i][j][1] = np.log10(xs[i][j]+xerrs[i][j][1])-np.log10(xs[i][j])
                    xs[i][j] = np.log10(xs[i][j])
            for j,ylist in enumerate(yarys[i]):
                ys[i][j] = self.avgfn(ylist)
                yerrs[i][j] = self.errfn(ylist)
                if (xs[i][j]>0 and ys[i][j]>0) or (logx and ys[i][j]>0):
                    empty = False
                if logy:
                    yerrs[i][j][0] = np.log10(ys[i][j])-np.log10(ys[i][j]-yerrs[i][j][0])
                    yerrs[i][j][1] = np.log10(ys[i][j]+yerrs[i][j][1])-np.log10(ys[i][j])
                    ys[i][j] = np.log10(ys[i][j])
        return xs,ys,xerrs,yerrs,empty

    def process_datapoints(self,plot_type,ion,xVar,xarys,yarys,average='default',logx='guess',logy='guess',\
                           visibility_threshold=None,**kwargs):
        logx,logy = self.should_take_logs_xy(ion,xVar,logx,logy,**kwargs)
        if average!='default':
            oldplots = self.plots
            self.setPlots(average,**kwargs)
        if average == 'covering_fraction':
            for i in range(len(yarys)):
                assert visibility_threshold is not None and plot_type in [0,1,2,3]
                yarys[i] = (yarys[i]>visibility_threshold).astype(int)
        if average == 'scatter':
            xs,ys,empty = self.process_scatter_points(xarys,yarys,logx,logy,**kwargs)
            xerrs,yerrs = None,None
        elif plot_type in [0,1,2,3]:
            xs,ys,yerrs,empty = self.process_errbars_onlyvertical(xarys,yarys,logx,logy,**kwargs)
            for i in range(len(xs)):
                yerrs[i] = np.transpose(yerrs[i])
            xerrs = None
        elif plot_type in [4]:
            xs,ys,xerrs,yerrs,empty = self.process_errbars_vertandhoriz(xarys,yarys,logx,logy,**kwargs)
            for i in range(len(xs)):
                yerrs[i] = np.transpose(yerrs[i])
                xerrs[i] = np.transpose(xerrs[i])
        if average!='default':
            self.setPlots(oldplots,**kwargs)
        return xs,ys,xerrs,yerrs,empty

    def get_title_and_axislabels(self,plot_type,ion,ion_name,xVar,xVar_name,replacement_title=None,extra_title='',\
                                 average='default',logx='guess',logy='guess',**kwargs):
        logx,logy = self.should_take_logs_xy(ion,xVar,logx,logy,**kwargs)
        if average == 'default':
            average = self.plots
        if replacement_title:
            title = replacement_title
        else:
            if plot_type == 0:
                title = "CGM Sightline Data"
            elif plot_type == 1:
                title = ion_name+" CGM Sightline Data"
            elif plot_type == 2:
                title = ion_name+" CGM Parameter Dependence"
            elif plot_type == 3:
                title = "Galaxy 2 Parameter Dependence"
            elif plot_type == 4:
                title = "Sightline 2 Variable Correlation"
        title+=" ("+str(average)+") "+extra_title

        if plot_type == 0 and ":" in ion[0] and ion[0].split(":")[1] != "cdens":
            ylabel = "fraction of ion in state"
        elif plot_type == 0 and (":" not in ion[0] or ion.split(":")[1] == "cdens"):
            ylabel = "col dens"
        elif ion in intensives:
            ylabel = intensiveslabels[ion]
        elif ion in param_xVars:
            ylabel = param_unit_labels[ion]
        elif self.plots == 'covering_fraction':
            ylabel = ion_name + ' covering fraction'
        else:
            ylabel = ion_name
        if logy and self.plots != 'covering_fraction':
            ylabel = 'log '+ylabel

        if xVar in intensives:
            xlabel = intensiveslabels[xVar]
        elif xVar in param_xVars:
            xlabel = param_unit_labels[xVar]
        elif xVar in sightline_xVars:
            xlabel = sightline_unit_labels[xVar]
        else:
            xlabel = xVar_name
        if logx:
            xlabel = 'log '+xlabel

        return xlabel,ylabel,title

    def plot_on_ax(self,ax,plot_type,xs,ys,xerrs,yerrs,xlabel,ylabel,title,labels=None,
                   average='default',dots=False,grid=False,linestyle='',
                   fmt=None,coloration=None,xlims='default',ylims='default',markersize='default',alpha = 1,
                   **kwargs):
        coloration = coloration or [None]*len(xs)
        future_colors = []
        if xerrs is None:
            xerrs=[None]*len(xs)
        if yerrs is None:
            yerrs=[None]*len(xs)
        if dots or plot_type in [3]:
            for i in range(len(xs)):
                fmt=fmt or 'o'
                if markersize=='default':
                    markersize=6
                color_store = ax.plot(xs[i],ys[i],marker = fmt,linestyle=linestyle,label=labels[i],
                                      color=coloration[i],markersize=markersize,alpha = alpha)
                future_colors.append(color_store[0].get_color())
        elif average == 'scatter':
            if markersize=='default':
                markersize=6
            for i in range(len(xs)):
                color_store = ax.plot(xs[i],ys[i],'o',label=labels[i],color=coloration[i],markersize=2,alpha = alpha)
                future_colors.append(color_store[0].get_color())
        else:
            fmtdict = {"mean":'.',"median_std":',',"covering_fraction":',',"stddev":'.',"median":"."}
            if not fmt:
                fmt = fmtdict[self.plots]
            for i in range(len(xs)):
                color_store = ax.errorbar(xs[i],ys[i],xerr=xerrs[i],yerr=yerrs[i],label=labels[i],ls=linestyle,\
                             color = coloration[i],fmt = fmt,capsize = 3,alpha = alpha)
                future_colors.append(color_store[0].get_color())
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if xlims != 'default':
            ax.set_xlim(xlims)        
        if ylims != 'default':
            ax.set_ylim(ylims)
        if grid:
            ax.grid()
        ax.legend()
        return future_colors

    def get_observational_data(self,include_observations,plot_type,ion,xVar,logx='guess',logy='guess',
                               rlims = None,lq=None,observationArray=None,**kwargs):
        logx,logy = self.should_take_logs_xy(ion,xVar,logx,logy,**kwargs)
        if include_observations is None:
            return None
        if lq is not None:
            labels = lq[0]
            observationArray = lq[3]
        if observationArray is None:
            observationArray = [self.currentObservationArray]
        if rlims is None:
            rlims = np.array([0.1,np.inf])
        elif rlims == "all":
            rlims = [0.0,np.inf]
        if xVar == 'rdivR':
            pass
            #self.constrain_current_Quasar_Array("Rvir_is_real",bins=['True'],changeArrayName=False)
        xerr_arys = None
        if plot_type==0:
            xarys,yarys = self.get_xy_type0(xVar,ion,rlims)
            _,yerr_arys = self.get_xy_type0(xVar,ion+':eb',rlims)
        elif plot_type==1:
            xarys,yarys = self.get_xy_type1(xVar,ion,observationArray,rlims)
            _,yerr_arys = self.get_xy_type1(xVar,ion+':eb',observationArray,rlims)
        elif plot_type==2:
            xarys,yarys = self.get_xy_type2(xVar,ion,observationArray,rlims)
            _,yerr_arys = self.get_xy_type2(xVar,ion+':eb',observationArray,rlims)
        elif plot_type==3:
            xarys,yarys = self.get_xy_type3(xVar,ion,observationArray)
            _,yerr_arys = self.get_xy_type3(xVar,ion+':eb',observationArray,rlims)
        elif plot_type==4:
            xarys,yarys = self.get_xy_type4(xVar,ion,observationArray,rlims)
            _,yerr_arys = self.get_xy_type4(xVar,ion+':eb',observationArray,rlims)
        xarys_detections = np.copy(xarys)
        yarys_detections = np.copy(xarys)
        yerr_arys_detections = np.copy(xarys)
        xarys_nondetections = np.copy(xarys)
        yarys_nondetections = np.copy(xarys)
        yerr_arys_nondetections = np.copy(xarys)
        
        for i in range(len(xarys)):
            if logy:
                for j in range(len(xarys[i])):
                    yarys[i][j] = np.log10(yarys[i][j])
                    yerr_arys[i][j] = np.log10(yerr_arys[i][j])
            if logx:
                for j in range(len(xarys[i])):
                    xarys[i][j] = np.log10(xarys[i][j])
            xarys_detections[i] = xarys[i][yerr_arys[i]>0]
            yarys_detections[i] = yarys[i][yerr_arys[i]>0]
            yerr_arys_detections[i] = yerr_arys[i][yerr_arys[i]>0]
            xarys_nondetections[i] = xarys[i][yerr_arys[i]<0]
            yarys_nondetections[i] = yarys[i][yerr_arys[i]<0]
            yerr_arys_nondetections[i] = yerr_arys[i][yerr_arys[i]<0]
        return xarys_detections,yarys_detections,yerr_arys_detections,xarys_nondetections,yarys_nondetections,yerr_arys_nondetections
    
    def plot_observational_data(self,ax,plot_type,include_observations,obs_data,xlabel,ylabel,title,labels=None,
               grid=False,linestyle='',obs_coloration=None,
               fmt=None,coloration=None,xlims='default',ylims='default',markersize='default',
               **kwargs):
        xarys_detections,yarys_detections,yerr_arys_detections,xarys_nondetections,yarys_nondetections,yerr_arys_nondetections = obs_data
        coloration = coloration or obs_coloration
        if coloration is None:
            coloration = [None]*len(xarys_detections)
        if include_observations == 'only':
            ax.cla()
        detections = ax.errorbar([],[],
                                     xerr=None,yerr=None,label='Werk et al. 2013',ls=linestyle,
                                     fmt = 's',capsize = 3, mec = 'k',ecolor = 'k',mfc='w')
        for i in range(len(xarys_detections)):
            color_store = ax.errorbar(xarys_detections[i],yarys_detections[i],
                                     xerr=None,yerr=yerr_arys_detections[i],ls=linestyle,
                                     color = coloration[i],fmt = 's',capsize = 3,alpha = .5)
            nondetectioncolor = color_store[0].get_color()
            ax.errorbar(xarys_nondetections[i][yerr_arys_nondetections[i]==-2],
                        yarys_nondetections[i][yerr_arys_nondetections[i]==-2],
                        xerr=None,yerr=.15,uplims=True,ls=linestyle,
                        mec = nondetectioncolor,ecolor = nondetectioncolor,fmt = 's',capsize = 3,mfc='w',alpha = .5)
            ax.errorbar(xarys_nondetections[i][yerr_arys_nondetections[i]==-3],
                        yarys_nondetections[i][yerr_arys_nondetections[i]==-3],
                        xerr=None,yerr=.15,lolims=True,ls=linestyle,
                        mec = nondetectioncolor,ecolor = nondetectioncolor,fmt = 's',capsize = 3,mfc='w',alpha = .5)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        if xlims != 'default':
            ax.set_xlim(xlims)        
        if ylims != 'default':
            ax.set_ylim(ylims)
        if grid:
            ax.grid()
        ax.legend()
        
    def plot_err(self,ion,xVar='rdivR',ax=None,show=False,include_observations = None,**kwargs):
        if not ax:
            print("Current constraints (name): "+self.currentQuasarArrayName)
            ax = plt.subplots(1)[1]
        plot_type = self.decide_plot_type(ion,xVar,**kwargs)
        ion,ion_name,xVar,xVar_name,labels = self.get_labels_from_ion(plot_type,ion,xVar,**kwargs)
        kwargs['labels']=labels
        xarys,yarys = self.get_sightline_xy_vals(plot_type,ion,xVar,**kwargs)
        xs,ys,xerrs,yerrs,empty = self.process_datapoints(plot_type,ion,xVar,xarys,yarys,**kwargs)
        xlabel,ylabel,title = self.get_title_and_axislabels(plot_type,ion,ion_name,xVar,xVar_name,**kwargs)
        if not empty:
            obs_colors = self.plot_on_ax(ax,plot_type,xs,ys,xerrs,yerrs,xlabel,ylabel,title,**kwargs)
            if include_observations:
                obs_data = self.get_observational_data(include_observations,plot_type,ion,xVar,**kwargs)
                self.plot_observational_data(ax,plot_type,include_observations,obs_data,xlabel,ylabel,title,obs_coloration = obs_colors,**kwargs)
            if show:
                plt.show()
        else:
            print("No values detected!")

    def definecolorbar(self, bar_type = 'HotCustom',**kwargs):
        from matplotlib.colors import LinearSegmentedColormap
        f= 256.0
        if bar_type not in ('HotCustom','RainbowCustom','BlackandWhite'):
            raise Exception("Not a ColorMap. Please try another one.")
        if bar_type == 'HotCustom':
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

            custom = LinearSegmentedColormap('HotCustom', cdict)
        elif bar_type == 'RainbowCustom':
            bowdict = {'red': ((0.0, 1.0, 1.0),
                                 (1e-9, 0.0, 0.0),
                                 (0.3, 0.5, 0.5),
                                 (0.6, 0.7, 0.7),
                                 (0.9, 0.8, 0.8),
                                 (1.0, 0.0, 0.0)),
                         'green': ((0.0, 1.0, 1.0),
                                   (1e-9, 0.0, 0.0),
                                   (0.3, 0.8, 0.8),
                                   (0.6, 0.7, 0.7),
                                   (0.9, 0.0, 0.0),
                                   (1.0, 0.0, 0.0)),
                         'blue': ((0.0, 1.0, 1.0),
                                  (1e-9, 1.0, 1.0),
                                  (0.3, 1.0, 1.0),
                                  (0.6, 0.0, 0.0),
                                  (0.9, 0.0, 0.0), 
                                  (1.0, 0.0, 0.0))}
            custom = LinearSegmentedColormap('RainbowCustom', bowdict)
        elif bar_type == 'BlackandWhite':
            bwdict = {'red': ((0.0, 1.0, 1.0),
                                 (1e-9, 0.8, 0.8),
                                 (0.3, 0.6, 0.6),
                                 (0.6, 0.4, 0.4),
                                 (0.9, 0.2, 0.2),
                                 (1.0, 0.0, 0.0)),
                         'green': ((0.0, 1.0, 1.0),
                                 (1e-9, 0.8, 0.8),
                                 (0.3, 0.6, 0.6),
                                 (0.6, 0.4, 0.4),
                                 (0.9, 0.2, 0.2),
                                 (1.0, 0.0, 0.0)),
                         'blue': ((0.0, 1.0, 1.0),
                                 (1e-9, 0.8, 0.8),
                                 (0.3, 0.6, 0.6),
                                 (0.6, 0.4, 0.4),
                                 (0.9, 0.2, 0.2),
                                 (1.0, 0.0, 0.0))}
            custom = LinearSegmentedColormap('BlackandWhite', bwdict)
        return custom
    
    def process_xy_vals_hist(self,ion,xs,ys,xVar='rdivR',tolerance=1e-5,weights=True,logx='guess',logy='guess',**kwargs):
        logx,logy = self.should_take_logs_xy(ion,xVar,logx,logy,**kwargs)
        unique_xs = self.combine_xs(xs,tolerance)
        retxs = xs*0.0
        for x_value in unique_xs:
            mask = self.values_within_tolerance(xs,x_value,tolerance,symmetric=False).astype(float)
            retxs+= (mask*x_value)
        if weights:
            weight = xs*0.0
            for i,x in enumerate(xs):
                weight[i] = 1/np.sum(self.values_within_tolerance(xs,x,tolerance,symmetric=False).astype(float))
            cbarlabel = "Fraction of lines for fixed %s"%(xVar)
        else:
            weight = xs*0.0+1.0
            cbarlabel = "Total number of lines"
        xs=retxs
        if logx:
            ys = ys[xs>0]
            xs = xs[xs>0]
            xs = np.log10(xs)
        if logy:
            xs = xs[ys>0]
            ys = ys[ys>0]
            ys = np.log10(ys)
        return xs,ys,weight,cbarlabel,len(xs)<=0
    
    def plot_on_ax_hist(self,ax,fig,xs,ys,xlabel,ylabel,title,weight,cbarlabel,ns = (42,15),\
                        xlims='default',ylims='default',**kwargs):
        hotcustom = self.definecolorbar(**kwargs)
        H, xedges, yedges = np.histogram2d(xs, ys, bins=ns,weights = weight)
        H = H.T
        X, Y = np.meshgrid(xedges, yedges)
        cms=ax.pcolormesh(X,Y, H, cmap=hotcustom)
        fig.colorbar(cms,label = cbarlabel,ax=ax)
        def get_new_axislim(current,divideBy = 200):
            x1 = current[0]
            x2 = current[1]
            return x1-(x2-x1)/divideBy,x2+(x2-x1)/divideBy
        if xlims=='default':
            ax.set_xlim(get_new_axislim(ax.get_xlim()))
        else:
            ax.set_xlim(xlims)
        if ylims=='default':
            ax.set_ylim(get_new_axislim(ax.get_ylim()))
        else:
            ax.set_ylim(ylims)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)


    #summary: plots histogram(s) of ion column density vs incrementing xvariable, 
    #         uses a color bar to show the percentage of sightlines for a certain column density at a specific x value
    #
    #    def plot_err(self, ion, quasarArray = None, xVar = "r", save_fig = False, \
    #             reset = False, labels = None,extra_title = "",rlims = None,\
    #             tolerance = 1e-5,dots = False,logx = False,average = None,logy = True
    def plot_hist(self,ion,xVar='rdivR',ax=None,fig=None,show=False,**kwargs):
        if isinstance(ion,tuple):
            ion_name = ion[1]
            ion = ion[0]
        elif isinstance(ion,str):
            ion_name=ion
        if isinstance(xVar,tuple):
            xVar_name = xVar[1]
            xVar = xVar[0]
        elif isinstance(ion,str):
            xVar_name=xVar
        if not ax:
            assert fig is None
            print("Current constraints (name): "+self.currentQuasarArrayName)
            fig,ax = plt.subplots(1)
        xs,ys = self.get_sightline_xy_vals(0,[ion],xVar,**kwargs)
        xs=xs[0]
        ys=ys[0]
        xs,ys,weight,cbarlabel,empty = self.process_xy_vals_hist(ion,xs,ys,**kwargs)
        xlabel,ylabel,title = self.get_title_and_axislabels(1,ion,ion_name,xVar,xVar_name,**kwargs)
        if not empty:
            self.plot_on_ax_hist(ax,fig,xs,ys,xlabel,ylabel,title,weight,cbarlabel,**kwargs)
            if show:
                plt.show()
        else:
            print("No values detected!")
            
    def sort_by_2D(self, criteria_x, criteria_y, bins_x = [0,np.inf], bins_y = [0,np.inf], \
                   atEnd_x = False, atEnd_y = False, splitEven_x = 0, splitEven_y = 0, \
                   reverse_x = False,reverse_y = False):
        
        #conduct a sort on one dimension and get the labels and bins for the other
        labels_y, bins_y, quasarbins_y, obs_bins_y = self.sort_by(criteria_y, 
                                                      bins_y,
                                                      atEnd = atEnd_y,
                                                      splitEven = splitEven_y,
                                                      reverse = ~reverse_y)
        labels_x, bins_x, _  ,_          = self.sort_by(criteria_x, 
                                                      bins_x, 
                                                      atEnd = atEnd_x, 
                                                      splitEven = splitEven_x,
                                                      reverse = False)
        if reverse_x:
            labels_x = reversearray(labels_x)
        quasarBins = np.zeros((len(labels_y),len(labels_x)), dtype = object)
        obsBins = np.zeros((len(labels_y),len(labels_x)), dtype = object)

        #sort the other dimension and store sorted quasarbins
        for i,qlist_y in enumerate(quasarbins_y):
            sorter = MultiSphereSorter(qlist_y)
            obs_sorter = MultiSphereSorter(obs_bins_y[i])
            _,_, quasarbins_x = sorter.sort(criteria_x,bins_x,atEnd = atEnd_x)
            _,_, obs_bins_x = obs_sorter.sort(criteria_x,bins_x,atEnd = False)
            if reverse_x:
                quasarbins_x = reversearray(quasarbins_x)
                obs_bins_x = reversearray(obs_bins_x)
            for j,qlist_x in enumerate(quasarbins_x):
                quasarBins[i][j] = qlist_x
                obsBins[i][j] = obs_bins_x[j]
        return labels_x,labels_y,bins_x,bins_y,quasarBins,obsBins                      
    
    def faberplot(self,plot_type,yVar,labels_x=None,labels_y=None,quasarArray=None,obsArray=None,lq2=None,criteria_legend=None,\
                  bins_legend=None,sharex=True,sharey=True,figsize='guess', save_fig = False,dpi=300, **kwargs):
        #after using sort_by_2d to get a set of labels and a 2d array of quasarspheres,
        #ask plot_err or plot_hist for completed plots of type given, for 
        #quasars in that cell of quasarArray, put them in subplots of a n by m subplots object
        #and show that plot
        print("Current constraints (name): "+self.currentQuasarArrayName)
        if lq2:
            labels_x,labels_y,_,_,quasarArray,obsArray = lq2
        else:
            assert labels_x is not None and labels_y is not None and quasarArray is not None
        rows = len(quasarArray)
        cols = len(quasarArray[0])
        if figsize=='guess':
            figsize=(min(15,5*cols),3.5*rows)
        xsize = figsize[0]
        ysize = figsize[1]
        fig, axes = plt.subplots(rows,cols,figsize = figsize,sharex=sharex,sharey=sharey)
        if rows == 1:
            axes = [axes]
        if cols == 1:
            newaxes = []
            for ax in axes:
                newaxes.append([ax])
            axes=newaxes
        oldQuasarArray = self.currentQuasarArray
        oldObsArray = self.currentObservationArray
        axlolims = np.array([np.inf,np.inf])
        axuplims = np.array([-np.inf,-np.inf])
        for r in range (0,rows):
            for c in range (0,cols):
                if 1:
                    self.currentQuasarArray = quasarArray[r][c]
                    self.currentObservationArray = obsArray[r][c]
                    if criteria_legend and len(self.currentQuasarArray)>0:
                        lq = self.sort_by(criteria_legend,bins_legend,**kwargs)
                    else:
                        lq=None
                    if plot_type=='err':
                        self.plot_err(yVar, ax = axes[r][c], show = False,lq=lq, **kwargs)
                    elif plot_type=='hist':
                        self.plot_hist(yVar, ax = axes[r][c], show = False,fig=fig, **kwargs)
                        
                        current_axlolims = np.array([axes[r][c].get_xlim()[0],axes[r][c].get_ylim()[0]])
                        current_axuplims = np.array([axes[r][c].get_xlim()[1],axes[r][c].get_ylim()[1]])
                        axlolims = np.amin((axlolims,current_axlolims),axis=0)
                        axuplims = np.amax((axuplims,current_axuplims),axis=0)
                        axes[r][c].set_xlim((axlolims[0],axuplims[0]))
                        axes[r][c].set_ylim((axlolims[1],axuplims[1]))
                    if r == 0:
                        axes[r][c].set_xlabel('')
                        axes[r][c].set_title(labels_x[c])
                    if r == rows-1:
                        axes[r][c].set_title('')
                    if r!=0 and r!= rows-1:
                        axes[r][c].set_xlabel('')
                        axes[r][c].set_title('')
                    if c == cols-1:
                        axes[r][c].set_ylabel('')
                        right_ax = axes[r][c].twinx()
                        right_ax.set_ylabel(labels_y[r])
                        right_ax.yaxis.set_ticks_position('none') 
                        right_ax.yaxis.set_ticklabels('') 
                    if c!=0 and c!=cols-1:
                        axes[r][c].set_ylabel('')

                #except Exception as e:
                else:
                    print(e)
        if not save_fig == False:
            fname = self.save_fig_filename(save_fig, yVar, **kwargs)
            plt.savefig(fname,dpi=dpi)
        self.currentQuasarArray = oldQuasarArray
        self.currentObservationArray = oldObsArray

    def save_fig_filename(self, save_fig, yVar, xVar = 'rdivR', **kwargs):
        if isinstance(save_fig, str):
            return 'quasarscan/plots/' + str + '.png'
        elif save_fig == True:
            filename = 'quasarscan/plots/' + yVar.replace('/','div') + '_vs_' + xVar.replace('/','div') + '_' + str(datetime.datetime.now())[:19] + '.png'
            filename = filename.replace(' ','')
            return filename
        else:
            print('Please use a valid file name or pass True.')
    
