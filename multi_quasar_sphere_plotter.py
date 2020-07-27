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

#summary: search directory for textfiles
#
#inputs: inquasarscan: if True, only look down one level. If false, look two.
#        loadonly: if not 'all' only load certain simulations (e.g. 'VELA')
#
#outputs: textfiles: list of names of textfiles
def get_all_textfiles(inquasarscan,loadonly = 'all'):
    
    #pathname by default starts in output
    if inquasarscan:
        path = "output"
        print('still true lol')
        dirs = os.listdir(path)
    else:
        path = "quasarscan/output"
        dirs = os.listdir(path)

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
    for folderName in dirs:
        if not folderName.startswith(".") and one_is_in_name(folderName,loadonly):
            folderPath = path + "/" + folderName
            folderDirs = os.listdir(folderPath)
            for fileName in folderDirs:
                if not fileName.startswith("."):
                    textfiles.append(os.path.join(folderPath,fileName))
    return textfiles

#summary: search directory for textfiles
#
#inputs: inquasarscan: if True, only look down one level. If false, look two.
#        loadonly: if not 'all' only load certain observations (e.g. 'CASBAH')
#
#outputs: textfiles: list of names of textfiles
def get_all_observations(inquasarscan,loadonly = 'all'):
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

#summary: put each list in a list of lists of ions in alphabetical/roman numeral order
#
#inputs: ions: a list of lists of ions
#        flat: if true, return only a single list of all the ions
#
#outputs: all the lists sorted, or sorted and reduced if 'flat'
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

#summary: reverse an array
#
#inputs: ary: an array
#
#outputs: the reversed array
def reversearray(ary):
    ary = list(ary)
    ary.reverse()
    return np.array(ary)

#these are a number of global lists and dictionaries which are checked against in various places
stringcriteria = ["ions","name","simname","version","code","simnum","Rvir_is_real","compaction_stage"]
intensives = ["Z","T","rho"]
intensiveslabels = {"Z":"avg metallicity","T":"avg temperature","rho":"avg density"}
intensivespositions = {"Z":-1,"T":-2,"rho":-3}
sightline_xVars = ["r","rdivR","rMpc","theta","phi"]
param_xVars = ["redshift","a0","Mvir","gas_Rvir","star_Rvir","dm_Rvir","sfr","ssfr","L_mag","Mstar","Mgas","Rvir"]
sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
           "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)","phi" \
           :"azimuthal viewing angle (rad)","rMpc":"r (Mpc)"}
param_unit_labels = {"redshift":"z","a0":"a","Rvir":'Virial radius (kpc)',"Mvir":"Virial Mass (Msun)",\
                    "gas_Rvir":"Gas Mass within Rvir (Msun)","Mgas":"Gas Mass within Rvir (Msun)","star_Rvir":"Stellar Mass within Rvir (Msun)",\
                    "Mstar":"Stellar Mass within Rvir (Msun)","dm_Rvir":"Dark Matter Mass within Rvir (Msun)","sfr":"Star Formation Rate (Msun yr-1)",\
                    "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)","L_mag":"Magnitude of Angular Momentum"}

class MultiQuasarSpherePlotter():
    #summary: initialize mq and load all data
    #
    #inputs: loadonly: if not 'all' only load certain simulations (e.g. 'VELA')
    #        loadobs: if not 'all' only load certain simulations (e.g. 'COS-Halos') [I'm not sure this works]
    #        textfiles: usually None. Could give a specific list of textfiles to use if you don't want to search here
    #        cleanup: delete textfiles that fail safety check (will ask user permission first)
    #        plots: default plots value to use. Will default to mean if not given
    #        throwErrors: if true, throw errors when reading textfiles if broken. If false, skip ones that create errors
    #        safetycheck: if False, skip safetychecks and just use whatever you load
    #
    #outputs: MultiQuasarSpherePlotter object, usually called 'mq'
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

    #summary: print length of list of QuasarSpheres
    #
    #inputs: None
    #
    #outputs: length of currentQuasarArray
    def length(self):
        return len(self.currentQuasarArray)
    
    #summary: print full list of QuasarSpheres, with metadata if requested (by 'criteria')
    #
    #inputs: *criteria: any criteria for a quasarSphere
    #
    #outputs: None, prints list of current quasarspheres
    def list_all_QuasarSpheres(self, *criteria,log=False):
        s = ""
        if len(self.currentQuasarArray)==0:
            print("no QuasarSpheres")
        for q in sorted(self.currentQuasarArray, key=lambda x: (x.name,x.rounded_redshift)):
            s += q.name +" at z=%s ("%q.rounded_redshift
            for c in criteria:
                v = eval('q.%s'%c)
                if isinstance(v,str) or isinstance(v,list):
                    s+="%s = %s, "%(c,v)
                elif log:
                    s+="%s = 10^%0.3f, "%(c,np.log10(v))
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
            
    #summary: tests each condition, each condition is a safety check on loaded data
    #
    #inputs: q: a quasarArray object
    #        require_metadata: if true, check that metadata loads accurately
    #
    #outputs: True if no problem, False otherwise
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

    #summary: cancel all constraints
    #
    #inputs: None
    #        
    #outputs: None, changes state of mq.currentQuasarArray and mq.currentQuasarArrayName
    def reset_current_Quasar_Array(self):
        self.currentQuasarArray = np.copy(self.quasarArray)
        self.currentObservationArray = np.copy(self.observationArray)
        self.currentQuasarArrayName = ''
        
    #summary: Check if you can constrain by a parameter
    #
    #inputs: constrainCriteria:  what criteria to constrain by
    #        atEnd: list of lists of quasarSphere objects that fit in the bins
    #        
    #outputs: oldQuasarArray: copy of currentQuasarArray that can be modified without affecting original
    #         if constraint is illegal, raises error
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

    #summary: helper for 'sort_by'. Sets up bins if needed
    #
    #inputs: criteria: list of strings for labelling points in legend
    #        bins: the bins to compare to 
    #        atEnd: list of lists of quasarSphere objects that fit in the bins
    #        
    #outputs: bins: the bins to compare to 
    #         oldQuasarArray: copy of currentQuasarArray that can be modified without affecting original
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
    
    #summary: helper for 'sort_by'. cleans up and gives options for bins 
    #
    #inputs: labels: list of strings for labelling points in legend
    #        bins: the bins to compare to 
    #        quasarBins: list of lists of quasarSphere objects that fit in the bins
    #        onlyNonempty: if True, remove (and do not plot) empty bins
    #        reverse: if True, return bins and quasarBins sorted in reverse order
    #        
    #outputs: labels: list of strings for labelling points in legend
    #         bins: the bins to compare to 
    #         quasarBins: list of lists of quasarSphere objects that fit in the bins
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
    
    
    #summary: splits currentQuasarArray into particular bins, either calculated on the fly or given
    #
    #inputs: criteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: a list of n numbers, which determine n-1 bins in between them. If string param, n strings which each
    #              constitute a bin
    #        atEnd: if True, compare values by their final (z=1 or z=minimum among remaining values) value, 
    #               not the current one
    #        splitEven: a number of bins to split into. The bins will be chosen so each has the same number of members.
    #        **kwargs: onlyNonempty,reverse ['postprocess_sorted']
    #        
    #outputs: labels: list of strings for labelling points in legend
    #         bins: the bins to compare to 
    #         quasarBins: list of lists of quasarSphere objects that fit in the bins
    #         obsBins: list of lists of observationalQuasarSphere objects that fit in the bins
    #          
    #         NOTE: These are usually combined together and considered an 'lq' object, passed directly 
    #            into most plots (any except type 0) general use case is e.g.
    #            >>>lq = mq.sort_by('Mvir',[0,10**11,np.inf])
    #            >>>mq.plot_err('O VI',lq=lq)
    def sort_by(self, criteria, bins = [0,np.inf],atEnd = False,splitEven = 0,**kwargs):
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

    #summary: similar to 'constrain_current_Quasar_Array' but specific to gasbins
    #
    #inputs: gasbintype: only keep lines where this type of gasbin was defined
    #        
    #outputs: none, changes state of 'currentQuasarArray'    
    def constrain_via_gasbins(self,gasbintype=None):
        if gasbintype == None:
            gasbintype = input("Available bins are: %s"%gasbinning.possible_bin_types)
        g = gasbinning.GasBinsHolder(bins=[gasbintype])
        toReturn = []
        for q in self.currentQuasarArray:
            if g.get_bin_str() in q.gasbins.get_bin_str():
                toReturn.append(q)
        self.currentQuasarArray = np.array(toReturn)
        
    #summary: helper function for constrain_array_helper
    #
    #inputs: constrainCriteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: either a list of two numbers, if a numerical criteria, or several strings if string criteria
    #              can leave as None if splitEven is used
    #        exclude: remove the ones listed, instead of keep the ones listed
    #        
    #outputs: bins: the bins to compare to 
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

    #summary: helper function for constrain_current_Quasar_Array
    #
    #inputs: sorter: MultiSphereSorter helper object, see 'sorter.py'
    #        obs_sorter: MultiSphereSorter helper object, see 'sorter.py'
    #        constrainCriteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: either a list of two numbers, if a numerical criteria, or several strings if string criteria
    #              can leave as None if splitEven is used
    #        splitEven: if 'high' or 'low' split into upper or lower half
    #        atEnd: if True, compare values by their final (z=1 or z=minimum among 
    #                 remaining values) value, not the current one
    #        set_main_array: make this array the natural basis which you will "reset" to (if working interactively)
    #        sortobs: if "default" only sort observations when constraining by 
    #                 things that are known for observations. Otherwise 
    #                 sort observations every time, filtering out for unknown data
    #        
    #outputs: bins: the bins to compare to   
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
    
    #summary: change the record of what constraints have been applied
    #
    #inputs: constrainCriteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: either a list of two numbers, if a numerical criteria, or several strings if string criteria
    #              can leave as None if splitEven is used
    #        changeArrayName: use False if you don't want to record a constraint (i.e. if they're getting too complicated)
    #        exclude: remove the ones listed, instead of keep the ones listed
    #        
    #outputs: None, changes state of mq.currentQuasarArrayName
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

    #summary: restricts to only quasarspheres with galaxy parameters within certain limits
    #
    #inputs: constrainCriteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: either a list of two numbers, if a numerical criteria, or several strings if string criteria
    #              can leave as None if splitEven is used
    #        **kwargs: changeArrayName, exclude ['change_array_name']
    #                  splitEven,atEnd,set_main_array,sortobs ['constrain_array_helper']
    #        
    #outputs: return the bins used (in the case of 'low' or 'high' for example it'll tell you the cutoff)
    def constrain_current_Quasar_Array(self, constrainCriteria,bins=None,**kwargs):
        self.check_criteria_works(constrainCriteria,**kwargs)
        bins = self.get_bin_values(constrainCriteria,bins,**kwargs)
        sorter = MultiSphereSorter(self.currentQuasarArray)
        obs_sorter = MultiSphereSorter(self.currentObservationArray)
        bins = self.constrain_array_helper(sorter,obs_sorter,constrainCriteria,bins,**kwargs)
        self.change_array_name(constrainCriteria,bins,**kwargs)
        return bins
    
    #summary: Sets mq's "plots" parameter. Several options, including median, mean, scatter, covering_fraction
    #         and combinations.
    #
    #inputs: plots: string choosing among ["mean","stderr","stddev","median","med","med_noquartiles",
    #                                      "median_std","scatter",'covering_fraction','cvf']
    #        quartiles: if plots == 'median' then you can give specific quartiles. Default is 40th % and 60th %
    #        
    #outputs: none. Change of state to mq 
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
            
    #summary: Fetches y values from saved data. First, has to get relevant info, then has 
    #         to perform math formula if any
    #
    #inputs: gq: GeneralizedQuasarSphere object (one or more quasarSpheres squished into one)
    #        stringVar: formula to show. Can use basic '+-*/()' characters. Variables have the names
    #                   'O VI' or 'O VI:temperature:hot' or things like this. Look at gasbinning.py for
    #                   more info
    #        
    #outputs: ys: unprocessed y data            
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

    #summary: In a type 0 plot, get sightline values for y and either discrete sightline values for x
    #         or metadata values for x. (it really just calls 'get_xy_type1' or 'get_xy_type2')
    #
    #inputs: xVar: xVar or xVar formula being displayed
    #        yVars: list of ions or formulas being displayed
    #        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
    #        
    #outputs: xs: raw x data for processing
    #         ys: raw y data for processing 
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

    #summary: In a type 1 plot, get sightline values for y and discrete sightline values for x
    #
    #inputs: xVar: xVar or xVar formula being displayed
    #        yVar: ion or formula being displayed
    #        quasarArray: array to use (from 'get_sightline_xy_vals')
    #        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
    #        
    #outputs: xs: raw x data for processing
    #         ys: raw y data for processing  
    def get_xy_type1(self,xVar,yVar,quasarArray,rlims):
        if rlims is None:
            rlims = [0.1,np.inf]
        elif rlims == "all":
            rlims = [0.0,np.inf]
        vardict = {"theta":1,"phi":2,"r":3,"rdivR":3,"rMpc":3}
        distances = "kpc" if xVar == "r" or xVar == "rMpc" else "Rvir"
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
            if xVar == "rMpc":
                xs[i]/=1000
            ys[i] = self.get_yVar_from_str(gq,yVar)
            rs = gq.info[:,3]
            acceptedLines = np.logical_and(rlims[0]<=rs,rs<=rlims[1])
            xs[i] = xs[i][acceptedLines]
            ys[i] = ys[i][acceptedLines]
        return xs,ys
    
    #summary: In a type 2 plot, get sightline values for y and metadata values for x
    #
    #inputs: xVar: xVar or xVar formula being displayed
    #        yVar: ion or formula being displayed
    #        quasarArray: array to use (from 'get_sightline_xy_vals')
    #        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
    #        
    #outputs: xs: raw x data for processing
    #         ys: raw y data for processing    
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
    
    #summary: In a type 3 plot, get values in the same way for both x and y. And they're just metadata table lookups
    #
    #inputs: xVar: xVar or xVar formula being displayed
    #        yVar: ion or formula being displayed
    #        quasarArray: array to use (from 'get_sightline_xy_vals')
    #        
    #outputs: xs: raw x data for processing
    #         ys: raw y data for processing
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
    
    #summary: In a type 4 plot, get values in the same way for both x and y
    #
    #inputs: xVar: xVar or xVar formula being displayed
    #        yVar: ion or formula being displayed
    #        quasarArray: array to use (from 'get_sightline_xy_vals')
    #        rlims: what impact parameters to accept. (from 'get_sightline_xy_vals')
    #        
    #outputs: xs: raw x data for processing
    #         ys: raw y data for processing
    def get_xy_type4(self, xVar, yVar, quasarArray, rlims):
        if rlims is None:
            rlims = [0.1,1.0]
        elif rlims == "all":
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
    
    #summary: helper for combine_xs. Checks if values are within tolerance
    #
    #inputs: x_in_list: value we're checking proximity to
    #        x_comp: value we're checking for proximity of
    #        tolerance: If x values are within tolerance, then count them as the same. 
    #        symmetric: if true, check values above and below. If false, only check above
    #                   (want to return from combine_xs the bottom end of the bins then)
    #        
    #outputs: True if close enough, False if not
    def values_within_tolerance(self,x_in_list,x_comp,tolerance,symmetric = True):
        tocompare = x_in_list-x_comp
        if symmetric:
            tocompare = np.abs(tocompare)
        return np.logical_and(-tolerance<=tocompare,tocompare<=tolerance)

    #summary: In a type 0, 1, or 2 plot, this is used to combine together a number of sightlines into 
    #         single data points affiliated with certain xVar values. This function only creates the 
    #         xVar values to use, they are applied in 'process_errbars_onlyvertical'
    #
    #inputs: x_variable: the set of all xVar numerical values
    #        tolerance: If x values are within tolerance, then count them as the same. 
    #        
    #outputs: unique_xs (for use in 'process_errbars_onlyvertical')
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
    
    #summary: given ion and xVar, figure out which kind of plot it is and 
    #         check the other information given is consistent with that type
    #
    #inputs: ion: ion or formula being displayed. Can be a tuple (formula, name)
    #        xVar: xVar or xVar formula being displayed. Can be a tuple (formula, name)
    #        labels: list of labels for all lines to plot. if labels, return them as is
    #        lq: this will include labels, if labels return as is
    #        
    #outputs: plot_type: 0,1,2,3,4 depending on kinds of variables
    #             - Type 0 compares several different y variables for the same set of sightlines. Each curve is a different 
    #               variable, and all variables are applied to the whole sample. Any xVar works.
    #             - Type 1 compares a single y variable for the same set of sightlines. Each curve is a different 
    #               subsample (i.e. a different galaxy snapshot/subset). Only sightline xVars accepted
    #             - Type 2 compares a single y variable for the same set of sightlines. Each curve is a different 
    #               subsample (i.e. a different galaxy subset). Only param xVars accepted (so probably only one datapoint per 
    #               galaxy)
    #             - Type 3 compares 2 params. Doesn't use our data really so not important. 
    #               Just lets us fit in Nir's data catalogue into our plotting system
    #             - Type 4 compares 2 ions or formulas against each other. Very useful for scatterplots.
    #               colored by subset or snapshot.
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

    #summary: In a type 0 plot, or a plot with no labels specified, guess what the labels should be
    #         also disaggregates ion formula and ion name if given
    #
    #inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        ion: ion or formula being displayed. Can be a tuple (formula, name)
    #        xVar: xVar or xVar formula being displayed. Can be a tuple (formula, name)
    #        labels: list of labels for all lines to plot. if labels, return them as is
    #        lq: this will include labels, if labels return as is
    #        
    #outputs: ion: formula to use to calculate y vals
    #         ion_name: name to use for those y vals
    #         xVar: formula to use to calculate x vals
    #         xVar_name: name to use for those y vals
    #         labels: list of list of labels for all lines to plot
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

    #summary: Get raw data for processing from saved data. Basically this just 
    #         calls one of four different functions depending on plot_type
    #
    #inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        ion: ion or formula being displayed
    #        xVar: xVar or xVar formula being displayed
    #        lq: tuple of splitting, from 'sort_by'
    #        quasarArray: array to use if not 'currentQuasarArray'
    #        rlims: what impact parameters to accept. Default is .1 Rvir to infinity
    #        
    #outputs: xarys: raw x data for processing
    #         yarys: raw y data for processing
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
            
        return xarys,yarys

    #summary: Detect whether x and y variables are naturally plotted as logs
    #
    #inputs: ion: ion or formula to use
    #        xVar: xVar or formula to use
    #        logx: whether to plot x logarithmically. If not 'guess', return whatever was given
    #        logy: whether to plot y logarithmically. If not 'guess', return whatever was given
    #        average: basically checking if you're plotting 'covering_fraction'
    #        
    #outputs: logx: whether to plot x logarithmically
    #         logy: whether to plot y logarithmically
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
    
    #summary: if not plotting errorbars but instead plotting random individual sightlines 
    #         data (average = 'scatter'), this will return the values needed
    #
    #inputs: xarys: raw x data (from 'get_sightline_xy_vals')
    #        yarys: raw y data (from 'get_sightline_xy_vals')
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        subsample: only plot this fraction of points. (half of points, 1/3 of points etc)
    #        offsetx: if True, calculate a small offset for x values and add it. (so errorbars don't overlap)
    #        tolerance: If x values are within tolerance, then count them as the same. The only effect on this is 
    #                   if offsetx is true
    #        
    #outputs: xs: processed x data to plot
    #         ys: processed y data to plot
    #         empty: if True, no data points to plot
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
    
    #summary: calculate the averages and errorbars in any plot except type 4
    #
    #inputs: xarys: raw x data (from 'get_sightline_xy_vals')
    #        yarys: raw y data (from 'get_sightline_xy_vals')
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        offsetx: if True, calculate a small offset for x values and add it. (so errorbars don't overlap)
    #        tolerance: If x values are within tolerance, then count them as the same. Might need to change for ssfr
    #        
    #outputs: xs: processed x data to plot
    #         ys: processed y data to plot
    #         yerrs: processed y error data to plot
    #         empty: if True, no data points to plot
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
            if isinstance(offsetx,tuple):
                xs[i]+=offsetx[1]
                offsetx_bool=offsetx[0]
            elif isinstance(offsetx,float) or isinstance(offsetx,int):
                xs[i]+=offsetx
                offsetx_bool = False
            elif isinstance(offsetx,bool):
                offsetx_bool = offsetx
            if offsetx_bool and len(xs[i])>0:
                scatterwidth = np.min(np.diff(np.unique(xs[i])))
                if logx:
                    xs[i][xs[i]>0]*=10**((float(i)/len(xs)-.5)*scatterwidth*.2)
                else:
                    xs[i]+=(float(i)/len(xs)-.5)*scatterwidth*.5

            if logx:
                xs[i] = np.log10(xs[i])
        return xs,ys,yerrs,empty
    
    #summary: calculate the averages and errorbars in type 4 plot
    #
    #inputs: xarys: raw x data (from 'get_sightline_xy_vals')
    #        yarys: raw y data (from 'get_sightline_xy_vals')
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        
    #outputs: xs: processed x data to plot
    #         ys: processed y data to plot
    #         xerrs: processed x error data to plot (if type 4, otherwise None)
    #         yerrs: processed y error data to plot
    #         empty: if True, no data points to plot
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

    #summary: helper function for plot_err that gets the appropriate data from the different data gathering functions
    #
    #inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        ion: ion or formula being displayed
    #        xVar: xVar or xVar formula being displayed
    #        xarys: raw x data (from 'get_sightline_xy_vals')
    #        yarys: raw y data (from 'get_sightline_xy_vals')
    #        average: what kind of average to use. See 'setPlots' for details. default is take mean and stderr
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        visibility_threshold: a float value. Count how many of the data points 
    #                              are above the value instead of averaging the values themselves. Only used
    #                              if average = 'covering_fraction'
    #        
    #outputs: xs: processed x data to plot
    #         ys: processed y data to plot
    #         xerrs: processed x error data to plot (if type 4, otherwise None)
    #         yerrs: processed y error data to plot
    #         empty: if True, no data points to plot
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
    
    #summary: helper function for plot_err that gets the labels of things
    #
    #inputs: plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        ion: ion or formula being displayed
    #        ion_name: display name for ion. Especially useful for formulas
    #        xVar: xVar or xVar formula being displayed
    #        xVar_name: display name for xVar. Especially useful for formulas and for Type 4 plots
    #        replacement_title: clear title and replace with this
    #        extra_title: add this to the autogenerated title
    #        average: what kind of average to use. See 'setPlots' for details. default is take mean and stderr
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        
    #outputs:  xlabel: label for x axis
    #          ylabel: label for y axis
    #          title: label for top of figure
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

    #summary: helper function for plot_err that handles simulation data
    #
    #inputs: ax: axis (from 'plot_err')
    #        plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        xs: x values (calculated in 'process_datapoints')
    #        ys: y values (calculated in 'process_datapoints')
    #        xerrs: x errors (calculated in 'process_datapoints'). Should be None unless type 4 plot
    #        yerrs: y values (calculated in 'process_datapoints')
    #        xlabel: label to write along x axis (from 'get_title_and_axislabels')
    #        ylabel: label to write along y axis (from 'get_title_and_axislabels')
    #        title: label to write along top (from 'get_title_and_axislabels')
    #        ylabel: label to write along y axis (from 'get_title_and_axislabels')
    #        labels: list of labels for curves
    #        average: what kind of average to use. See 'setPlots' for details. default is take mean and stderr
    #        dots: If True, plot points instead of errorbars
    #        grid: whether to plot grid underneath values
    #        linestyle: how to connect points, default is no connection
    #        linewidth: thickness of connecting line
    #        fmt: style of points ('o','s','x','v', etc)
    #        coloration: a list of colors. If none go through the default matplotlib colors
    #        xlims: limits of x axis of plot. Default is matplotlib default
    #        ylims: limits of y axis of plot. Default is matplotlib default
    #        markersize: size of datapoints on plot
    #        alpha: transparency of datapoints
    #        elinewidth: thickness of errorbars. Default is linewidth
    #        
    #outputs:  future_colors: list of colors plotted to replicate in other plots if desired
    def plot_on_ax(self,ax,plot_type,xs,ys,xerrs,yerrs,xlabel,ylabel,title,labels=None,
                   average='default',dots=False,grid=False,linestyle='',linewidth = 1.5,
                   fmt=None,coloration=None,xlims='default',ylims='default',markersize='default',
                   alpha = 1.0,elinewidth=None,
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
                                      color=coloration[i],markersize=markersize,alpha = alpha,linewidth=linewidth)
                future_colors.append(color_store[0].get_color())
        elif average == 'scatter':
            if markersize=='default':
                markersize=6
            for i in range(len(xs)):
                color_store = ax.plot(xs[i],ys[i],'o',label=labels[i],color=coloration[i],markersize=markersize,
                                        alpha = alpha,linewidth=linewidth)
                future_colors.append(color_store[0].get_color())
        else:
            fmtdict = {"mean":'.',"median_std":',',"covering_fraction":',',"stddev":'.',"median":"."}
            if not fmt:
                fmt = fmtdict[self.plots]
            if fmt=='fill':
                for i in range(len(xs)):
                        color_store = ax.plot(xs[i],ys[i],label=labels[i],ls=linestyle,\
                                     color = coloration[i],marker = fmtdict[self.plots],alpha = alpha, linewidth=linewidth)
                        future_colors.append(color_store[0].get_color())
                        ax.fill_between(xs[i],ys[i]-yerrs[i][0,:],ys[i]+yerrs[i][1,:],color = future_colors[-1],\
                                        alpha = alpha/2,linewidth = linewidth/2)
            else:
                for i in range(len(xs)):
                        color_store = ax.errorbar(xs[i],ys[i],xerr=xerrs[i],yerr=yerrs[i],label=labels[i],ls=linestyle,\
                                     color = coloration[i],fmt = fmt,capsize = 3,alpha = alpha, linewidth=linewidth, \
                                     elinewidth=elinewidth)
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

    #summary: gets processed data for 'plot_observational_data' from the current observationarray
    #
    #inputs: include_observations: if None return None
    #        plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        ion: the ion or formula to plot(can only use single ion b/c shows detailed ion info as color)
    #        xVar: the xVar to plot. Can use any variable, will figure out what kind in 'decide_plot_type'
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        rlims: what impact parameters to include. Default is .1 Rvir to infinity
    #        lq: tuple of splitting, from 'sort_by'
    #        observationArray: observations to use if not all current observations
    #
    #outputs:  None. This actually does the plot operation (ax.errorbar)
    def get_observational_data(self,include_observations,plot_type,ion,xVar,logx='guess',logy='guess',
                               rlims = None,lq=None,observationArray=None,**kwargs):
        def flatten_if_needed(ary):
            try:
                for item in ary:
                    len(item)
                return np.concatenate(ary)
            except:
                return ary
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
        xarys_nondetections_up = np.copy(xarys)
        xarys_nondetections_down = np.copy(xarys)
        yarys_nondetections_up = np.copy(xarys)
        yarys_nondetections_down = np.copy(xarys)
        for i in range(len(xarys)):
            xarys_detections[i] = xarys[i][yerr_arys[i]>=0]
            yarys_detections[i] = yarys[i][yerr_arys[i]>=0]
            yerr_arys_detections[i] = yerr_arys[i][yerr_arys[i]>=0]
            xarys_nondetections_down[i] = xarys[i][yerr_arys[i]==-2]
            yarys_nondetections_down[i] = yarys[i][yerr_arys[i]==-2]
            xarys_nondetections_up[i] = xarys[i][yerr_arys[i]==-3]
            yarys_nondetections_up[i] = yarys[i][yerr_arys[i]==-3]
            if logy:
                yarys_detections[i] = np.log10(flatten_if_needed(yarys_detections[i]))
                yerr_arys_detections[i] = np.log10(flatten_if_needed(yerr_arys_detections[i]))
                yarys_nondetections_down[i] = np.log10(flatten_if_needed(yarys_nondetections_down[i]))
                yarys_nondetections_up[i] = np.log10(flatten_if_needed(yarys_nondetections_up[i]))
            if logx:
                xarys_detections[i] = np.log10(flatten_if_needed(xarys_detections[i]))
                xarys_nondetections_down[i] = np.log10(flatten_if_needed(xarys_nondetections_down[i]))
                xarys_nondetections_up[i] = np.log10(flatten_if_needed(xarys_nondetections_up[i]))
        return xarys_detections,yarys_detections,yerr_arys_detections,xarys_nondetections_up,yarys_nondetections_up,xarys_nondetections_down,yarys_nondetections_down
    
    #summary: helper function for plot_err that handles observational data
    #
    #inputs: ion: the ion or formula to plot(can only use single ion b/c shows detailed ion info as color)
    #        xVar: the xVar to plot. Can use any variable, will figure out what kind in 'decide_plot_type'
    #        ax: axis (from 'plot_err')
    #        fig: fig (from 'plot_err')
    #        plot_type: 0,1,2,3,4 depending on kinds of variables (see 'decide_plot_type')
    #        include_observations: 'only' or 'both' (if False, this func won't be called)
    #        obs_data: package of data (from 'get_observational_data')
    #        xlabel: label to write along x axis (from 'get_title_and_axislabels')
    #        ylabel: label to write along y axis (from 'get_title_and_axislabels')
    #        title: label to write along top (from 'get_title_and_axislabels')
    #        labels: list of labels for curves
    #        grid: whether to plot grid underneath values
    #        linestyle: how to connect points, default is no connection (which is probably best in all cases for observations)
    #        obs_coloration: if observations should be colored different from simulation, this is a list of colors
    #        alpha: transparency of observations
    #        fmt: style of points ('o','s','x','v', etc)
    #        coloration: a list of colors, inherited from simulation colors. If obs_coloration, ignore this
    #        xlims: limits of x axis of plot. Default is matplotlib default
    #        ylims: limits of y axis of plot. efault is matplotlib default
    #
    #outputs:  None. This actually does the plot operation (ax.errorbar)
    def plot_observational_data(self,ax,plot_type,include_observations,obs_data,xlabel,ylabel,title,labels=None,
               grid=False,linestyle='',obs_coloration=None,alpha=.5,
               fmt=None,coloration=None,xlims='default',ylims='default',
               **kwargs):
        xarys_detections,yarys_detections,yerr_arys_detections,xarys_nondetections_up,yarys_nondetections_up,xarys_nondetections_down,yarys_nondetections_down = obs_data
        coloration = coloration or obs_coloration
        if coloration is None:
            coloration = [None]*len(xarys_detections)
        if include_observations == 'only':
            ax.cla()
        detections = ax.errorbar([],[],
                                     xerr=None,yerr=None,label='Werk 2013',ls=linestyle,
                                     fmt = 's',capsize = 3, mec = 'b',ecolor = 'b',mfc='b')
        for i in range(len(xarys_detections)):
            color_store = ax.errorbar(xarys_detections[i],yarys_detections[i],
                                     xerr=None,yerr=yerr_arys_detections[i],ls=linestyle,
                                     color = coloration[i],fmt = 's',capsize = 3,alpha = alpha)
            nondetectioncolor = color_store[0].get_color()
            ax.errorbar(xarys_nondetections_up[i],yarys_nondetections_up[i],xerr=None,yerr=.15,lolims=True,ls=linestyle,\
                         mec = nondetectioncolor,ecolor = nondetectioncolor,fmt = 's',capsize = 3,mfc='w',alpha = alpha)
            ax.errorbar(xarys_nondetections_down[i],yarys_nondetections_down[i],xerr=None,yerr=.15,uplims=True,ls=linestyle,\
                         mec = nondetectioncolor,ecolor = nondetectioncolor,fmt = 's',capsize = 3,mfc='w',alpha = alpha)
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
        
    #summary: plots errorbar(s) of ion column density vs incrementing xvariable, 
    #         uses a color bar to show the average (or median) of sightlines for a 
    #         certain column density at a specific x value
    #
    #inputs: ion: the ion or formula to plot(can only use single ion b/c shows detailed ion info as color)
    #        xVar: the xVar to plot. Can use any variable, will figure out what kind in 'decide_plot_type'
    #        ax: axis if exists
    #        fig: fig if exists
    #        show: True if show and clear fig
    #        include_observations: True if you want to plot observations alongside simulation data
    #        **kwargs: obs_data, obs_coloration,xlabel,ylabel,title,labels,
    #                  grid,linestyle,alpha,fmt,coloration,xlims,ylims ['plot_observational_data']
    #                  logx, logy, rlims, lq, observationArray ['get_observational_data']
    #                  average, dots, linewidth, markersize, elinewidth ['plot_on_ax']
    #                  replacement_title, extra_title ['get_title_and_axislabels']
    #                  visibility_threshold ['process_datapoints']
    #                  offsetx, tolerance ['process_errbars_onlyvertical']
    #outputs: None. Will plot errorbars if show is true, otherwise will add errorbar to plt
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
            
    #summary: a number of custom colormaps
    #
    #inputs: bar_type: one of three options: 'HotCustom','RainbowCustom','BlackandWhite'
    #        call 'plot_hist' with the keyword 'bar_type' to use any besides HotCustom
    #
    #outputs: colormap object
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
    
    #summary: creates data to plot for histogram
    #
    #inputs: ion: ion or formula to plot
    #        xs: raw x values for all sightlines(from 'get_sightline_xy_vals')
    #        ys: raw y values for all sightlines(from 'get_sightline_xy_vals')
    #        xVar: what the raw x values represent
    #        tolerance: how close values can be to count as a single value (can be off by floating point error)
    #                   in the case of small values (e.g. ssfr) this may need to be reduced
    #        weights: whether to calculate in terms of raw sightline numbers (False) or 
    #                 fraction of each unique xVar (True, default)
    #        weight: multiplier for values to scale by individual densities at each xVar 
    #                (calculated in 'process_xy_vals_hist')
    #        logx: whether to plot x logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #        logy: whether to plot y logarithmically. Default is "guess", it'll ask 'should_take_logs_xy'
    #
    #outputs: xs: x values for histogram
    #         ys: x values for histogram
    #         weight: weight values for each x in histogram
    #         cbarlabel: basically, tell whether or not values are weighted
    #         len(xs)<=0: true if the histogram will be empty (if so, it'll be skipped)
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

    #summary: helper function for plot_hist where plot function is actually called
    #
    #inputs: ax: matplotlib axis (from 'plot_hist')
    #        fig: matplotlib fig (from 'plot_hist')
    #        xs: x values (calculated in 'process_xy_vals_hist')
    #        ys: y values (calculated in 'process_xy_vals_hist')
    #        xlabel: label to write along x axis (from 'get_title_and_axislabels')
    #        ylabel: label to write along y axis (from 'get_title_and_axislabels')
    #        title: label to write top (from 'get_title_and_axislabels')
    #        weight: multiplier for values to scale by individual densities at each xVar 
    #                (calculated in 'process_xy_vals_hist')
    #        cbarlabel: label to write along colorbar (from 'get_title_and_axislabels')
    #        ns: grid spacing on plot. increase first number to get more, thinner boxes. Increase second number
    #            to get more, shorter boxes. Too many gives essentially binary results. Default usually works
    #        xlims: x limits, defaults to the max and min
    #        ylims: y limits, defaults to the max and min
    #
    #outputs: None. Actually does the plot operation (ax.pcolormesh)
    def plot_on_ax_hist(self,ax,fig,xs,ys,xlabel,ylabel,title,weight,cbarlabel,ns = (42,15),\
                        xlims='default',ylims='default',**kwargs):
        cmap = self.definecolorbar(**kwargs)
        H, xedges, yedges = np.histogram2d(xs, ys, bins=ns,weights = weight)
        H = H.T
        X, Y = np.meshgrid(xedges, yedges)
        cms=ax.pcolormesh(X,Y, H, cmap=cmap)
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
    #inputs: ion: the ion or formula to plot(can only use single ion b/c shows detailed ion info as color)
    #        xVar: the xVar to plot (can only use sightline params here)
    #        ax: axis if exists
    #        fig: fig if exists
    #        show: True if show and clear fig
    #        **kwargs: ns, xlims, ylims ['plot_on_ax_hist']
    #                tolerance, weights, logx, logy ['process_xy_vals_hist']
    #                bar_type ['definecolorbar']
    #
    #outputs: None. Will plot histogram if show is true, otherwise will add histogram to plt
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

    #summary: splits currentQuasarArray into two dimensions of bins, either calculated on the fly or given
    #
    #inputs: criteria_x: what criteria to constrain by horizontally (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        criteria_y: what criteria to constrain by vertically 
    #        bins_x: a list of n numbers, which determine n-1 horizontal bins in between them. If string  
    #              param, n strings which each constitute a bin
    #        bins_y: a list of n numbers, which determine n-1 horizontal bins in between them. If string  
    #              param, n strings which each constitute a bin
    #        atEnd_x: if True, compare horizontal values by their final (z=1 or 
    #                 z=minimum among remaining values) value, not the current one
    #        atEnd_y: if True, compare vertical values by their final (z=1 or 
    #                 z=minimum among remaining values) value, not the current one
    #        splitEven_x: a number of bins to split into horizontally. 
    #                     The bins will be chosen so each has the same number of members.
    #        splitEven_y: a number of bins to split into vertically. 
    #                     The bins will be chosen so each has the same number of members.
    #        reverse_x: if True, return x bins in reverse order
    #        reverse_y: if True, return y bins in reverse order
    #        
    #outputs: labels_x: list of strings for labelling panels on top of faberplot
    #         labels_y: list of strings for labelling panels on left of faberplot
    #         bins_x: the bins to compare to horizontally
    #         bins_y: the bins to compare to vertically
    #         quasarBins: 2D list of lists of quasarSphere objects that fit in the bins
    #         obsBins: 2D list of lists of observationalQuasarSphere objects that fit in the bins
    #          
    #         NOTE: These are usually combined together and considered an 'lq2' object, passed directly 
    #            into faberplots (any except type 0) general use case is e.g.
    #            >>>lq2 = mq.sort_by2D('Mvir','redshift',[0,10**11,np.inf],splitEven_y = 3)
    #            >>>mq.faberplot('err','O VI',lq=lq)       
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

    #summary: creates 2d array of labelled miniplots
    #
    #inputs: plot_type: either 'err' or 'hist', will call those above
    #        yVar: what to plot (also pass xVar, but it counts in **kwargs here)
    #        lq2: output of 'sort_by_2D'
    #              param, n strings which each constitute a bin
    #        criteria_legend: to make type 1,2,3,4 error plots, you need to pass in a way of sorting
    #                         the subpanels. So pass in a 'sort_by' criteria. But don't sort in advance
    #                         the sort will happen inside faberplot
    #        bins_legend: bins for sorting inside subpanels. [bins for legend, not legend for bins]
    #        sharex: if True, all subpanels share same x values
    #        sharey: if True, all subpanels share same y values
    #        figsize: size of finalized image. Will guess according to number of bins if none given
    #        save_fig: whether to save. if string, this will be the file's name
    #        dpi: resolution of saved fig
    #        
    #outputs: None, makes plot. Can save if 'save_fig'
    def faberplot(self,plot_type,yVar,lq2=None,criteria_legend=None,bins_legend=None,\
                  sharex=True,sharey=True,figsize='guess', save_fig = False,dpi=300,\
                  **kwargs):
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

    #summary: create name for saved fig if not given
    #
    #inputs: save_fig: either True or a string. if a string skip this func 
    #        yVar: y variable plotted
    #        xVar x variable plotted
    #        
    #outputs: returns a string of a filename
    def save_fig_filename(self, save_fig, yVar, xVar = 'rdivR', **kwargs):
        if isinstance(save_fig, str):
            return 'quasarscan/plots/' + save_fig + '.png'
        elif save_fig == True:
            filename = 'quasarscan/plots/' + yVar.replace('/','div') + '_vs_' + \
                        xVar.replace('/','div') + '_' + str(datetime.datetime.now())[:19] + '.png'
            filename = filename.replace(' ','')
            return filename
        else:
            print('Please use a valid file name or pass True.')
    
