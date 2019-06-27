import numpy as np
import os
import sys
import matplotlib.pyplot as plt
try:
    from quasarscan import quasar_sphere
    from quasarscan import ion_lists
    from quasarscan import gasbinning
    from quasarscan import roman
    level = 0
except:
    import quasar_sphere
    import ion_lists
    import gasbinning
    import roman
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

stringcriteria = ["ions","simname","version","code","simnum","Rvir_is_real"]
intensives = ["Z","T","rho"]
intensiveslabels = {"Z":"avg log metallicity","T":"avg log temperature","rho":"avg log density"}
intensivespositions = {"Z":-1,"T":-2,"rho":-3}
sightline_xVars = ["r","rdivR","theta","phi"]
param_xVars = ["redshift","a0","Mvir","gas_Rvir","star_Rvir","dm_Rvir","sfr","ssfr","L_mag","Mstar","Mgas"]
sightline_unit_labels = {"r":"r (kpc)","r>0":"r (kpc)","rdivR":"r/Rvir","rdivR>0":"r/Rvir",\
           "theta":"viewing angle (rad)","theta_r>0":"viewing angle (rad)","phi" \
           :"azimuthal viewing angle (rad)"}
param_unit_labels = {"redshift":"z","a0":"a","Mvir":"Virial Mass (Msun)",\
                    "gas_Rvir":"Gas Mass within Rvir (Msun)","Mgas":"Gas Mass within Rvir (Msun)","star_Rvir":"Stellar Mass within Rvir (Msun)",\
                    "Mstar":"Stellar Mass within Rvir (Msun)","dm_Rvir":"Dark Matter Mass within Rvir (Msun)","sfr":"Star Formation Rate (Msun yr-1)",\
                    "ssfr":"Specific Star Formation Rate (Msun yr-1 Mstar-1)","L_mag":"Magnitude of Angular Momentum"}

class MultiQuasarSpherePlotter():
    #param: textfiles     if a list of textfiles is specified, those specific textfiles will be loaded; else,
    #                     all textfiles in output are loaded
    def __init__(self, loadonly = "all",textfiles = None, cleanup = False,plots = "mean",throwErrors = False,safetycheck = True):
        self.plots = "mean"
        self.avgfn = np.mean
        self.setPlots(plots)
        self.quasarArray = []
        if textfiles is None:
            textfiles = get_all_textfiles(level,loadonly = loadonly)
        for textfile in textfiles:
            try:
                readvalsoutput = quasar_sphere.read_values(textfile)
                q = quasar_sphere.QuasarSphere(readvalsoutput = readvalsoutput)
                if safetycheck and self.pass_safety_check(q):
                    self.quasarArray.append(q)
                elif cleanup:
                    todo = raw_input("file %s did not pass safety check. Remove it? (y/n)"%textfile).lower()
                    os.remove(textfile) if todo == 'y' else None
            except Exception as e:
                print(textfile + " could not load because:")
                print(e)
                if throwErrors:
                    raise e
        self.quasarArray = np.array(self.quasarArray)
        self.observedQuasarArray = np.array(self.quasarArray)
        self.currentQuasarArray = np.copy(self.quasarArray)
        self.currentQuasarArrayName = ''
        
        if len(self.currentQuasarArray) == 0:
            print ("There are no quasarspheres stored in currentQuasarArray!")
    
    def length(self):
        return len(self.currentQuasarArray)
    def list_all_QuasarSpheres(self, *criteria):
        s = ""
        if len(self.currentQuasarArray)==0:
            print "no QuasarSpheres"
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
        print s[:-1]
            
    #tests each condition, each condition is a safety check
    def pass_safety_check(self, q,require_metadata = False):
        minlength = 10
        minions = []
        if q.length_reached < minlength:
            print "Length for %s at redshift %s is not valid." %(q.name, str(q.rounded_redshift))
            return False
        elif not all(x in q.ions for x in minions):
            print "Not all necessary ions present in %s at redshift %s."%(q.name, str(q.rounded_redshift))
            return False
        if require_metadata:
            if q.final_a0 is None:
                print "metadata for %s at redshift %s is not valid." %(q.name,str(q.rounded_redshift))
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
        if isinstance(bins, float) or isinstance(bins, int) or \
            (len(bins) == 1 and isinstance(bins[0], float)) or (len(bins) == 1 and isinstance(bins[0], int)):
            if isinstance(bins, list) or isinstance(bins,np.ndarray):
                bins = bins[0]
            bins = np.array([0.0, bins, np.inf])
        elif isinstance(bins, str) and criteria in stringcriteria:
             bins = [bins]
        if isinstance(atEnd,float):
            self.constrain_current_Quasar_Array("final_a0",[atEnd,np.inf],changeArrayName=False)
            if len(self.currentQuasarArray) == 0:
                print "No galaxies get to that high of a0"
            atEnd = True
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
    
    def constrain_current_Quasar_Array(self, constrainCriteria, bins=None, exploration_mode = False,atEnd = False,\
        splitEven = None,changeArrayName = True, set_main_array = False, exclude = False):
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
        if exclude:
            possible_bins = []
            for q in self.currentQuasarArray:
                possible_bins.append(eval("q.%s"%constrainCriteria))
            possible_bins = set(possible_bins)
            excluded_bins = bins
            bins = list(possible_bins.difference(set(bins)))
        if isinstance(atEnd,float):
            self.constrain_current_Quasar_Array("final_a0",[atEnd,np.inf],changeArrayName=False)
            if len(self.currentQuasarArray) == 0:
                print "No galaxies get to that high of a0"
            atEnd = True
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
        if set_main_array:
            self.quasarArray = np.copy(self.currentQuasarArray)
        
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
                    bins = excluded_bins
                for acceptedValue in bins:
                    self.currentQuasarArrayName += acceptedValue.replace(" ","")
        return bins

    def setPlots(self,plots,quartiles = None):
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
        def bylength(word1,word2):
            return len(word2)-len(word1)
        def sortlist(a):
            a.sort(cmp=bylength)
            return a
        strings_to_find = split_by_ops(stringVar)
        new_str_to_eval = stringVar
        strings_to_replace_with = {}
        for i,s in enumerate(strings_to_find):
            if len(s)>0:
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

    def values_within_tolerance(self,x_in_list,x_comp,tolerance,symmetric = True):
        tocompare = x_in_list-x_comp
        if symmetric:
            tocompare = np.abs(tocompare)
        return np.logical_and(-1e-5<=tocompare,tocompare<=tolerance)

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
            x_values.append(np.min([currentList]))
        if i == len(x_values_not_averaged)-1:
            x_values.append(x_values_not_averaged[-1]) 
        return np.array(x_values)

    def decide_plot_type(self,ion,xVar='rdivR',labels=None,quasarArray=None,lq=None,**kwargs):
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

    def get_labels_from_ion(self,plot_type,ion,labels=None,lq=None,**kwargs):
        if lq is not None:
            labels=lq[0]
        if isinstance(ion,tuple):
            ion_name=ion[1]
            ion = ion[0]
        else:
            ion_name=ion
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
        return ion,ion_name,labels

    def get_sightline_xy_vals(self,plot_type,ion,xVar = 'rdivR',lq=None,quasarArray=None,rlims=None,**kwargs):
        if lq is not None:
            quasarArray = lq[2]
        if quasarArray is None:
            quasarArray = [self.currentQuasarArray]
        if rlims is None:
            rlims = np.array([0.1,np.inf])
        elif rlims == "all":
            rlims = [0.0,np.inf]
        
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
        return xarys,yarys
    
    def should_take_logs_xy(self,ion,xVar,logx,logy,average='default',**kwargs):
        probablylinear = ['redshift','rounded_redshift','a0','Rvir']
        if isinstance(xVar,tuple):
            xVar = xVar[0]
        if logx=='guess':
            if xVar in sightline_xVars:
                logx=False
            elif xVar in param_xVars and xVar not in probablylinear:
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
            if offsetx:
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
                    yerrs[i][j][1] = np.log10(ys[i][j]+yerrs[i][j][1])-np.log10(ys[i][j])
                    ys[i][j] = np.log10(ys[i][j])
            xs[i] = xs[i][mask]
            ys[i] = ys[i][mask]
            yerrs[i] = yerrs[i][mask]
            if offsetx:
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
                if xs[i][j]>0 and ys[i][j]>0:
                    empty = False
                if logy:
                    yerrs[i][j][0] = np.log10(ys[i][j])-np.log10(ys[i][j]-yerrs[i][j][0])
                    yerrs[i][j][1] = np.log10(ys[i][j]+yerrs[i][j][1])-np.log10(ys[i][j])
                    ys[i][j] = np.log10(ys[i][j])
        return xs,ys,xerrs,yerrs,empty

    def process_datapoints(self,plot_type,ion,xarys,yarys,xVar='rdivR',average='default',logx='guess',logy='guess',\
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

    def get_title_and_axislabels(self,plot_type,ion,ion_name,xVar='rdivR',replacement_title=None,extra_title='',\
                                 average='default',logx='guess',logy='guess',**kwargs):
        logx,logy = self.should_take_logs_xy(ion,xVar,logx,logy,**kwargs)
        if isinstance(xVar,tuple):
            xVar = xVar[0]
            xVar_name = xVar[1]
        else:
            xVar_name = xVar
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
        title+=" ("+average+") "+extra_title

        if plot_type == 0 and ":" in ion[0] and ion[0].split(":")[1] != "cdens":
            ylabel = "fraction of ion in state"
        elif plot_type == 0 and (":" not in ion[0] or ion.split(":")[1] == "cdens"):
            ylabel = "col dens"
        elif ion in intensives:
            ylabel = intensiveslabels[ion]
        elif ion in param_xVars:
            ylabel = param_unit_labels[ion]
        elif ":" in ion and ion.split(":")[1] != "cdens":
            ylabel = "fraction of ion in state: %s"%(ion_name)
        elif ":" not in ion or ion.split(":")[1] == "cdens":
            ylabel = "col dens %s"%(ion_name)
        elif self.plots == 'covering_fraction':
            ylabel = ion_name + ' covering fraction'
        if logy and self.plots != 'covering_fraction':
            ylabel = 'log '+ylabel

        if xVar in intensives:
            xlabel = intensiveslabels[xVar]
        elif xVar in param_xVars:
            xlabel = param_unit_labels[xVar]
        elif xVar in sightline_xVars:
            xlabel = sightline_unit_labels[xVar]
        elif ":" in xVar and xVar.split(":")[1] != "cdens":
            xlabel = "fraction of xVar in state: %s"%(xVar_name)
        elif ":" not in xVar or xVar.split(":")[1] == "cdens":
            xlabel = "col dens %s"%(xVar_name)
        if logx:
            xlabel = 'log '+xlabel

        return xlabel,ylabel,title

    def plot_on_ax(self,ax,plot_type,xs,ys,labels,xerrs,yerrs,xlabel,ylabel,title,average='default',dots=False,fmt=None,coloration=None,**kwargs):
        coloration = coloration or [None]*len(xs)
        if xerrs is None:
            xerrs=[None]*len(xs)
        if yerrs is None:
            yerrs=[None]*len(xs)
        fmtdict = {"mean":'.',"median_std":',',"covering_fraction":',',"stddev":'.',"median":"."}
        if not fmt:
            fmt = fmtdict[self.plots]
        if dots or plot_type in [3]:
            for i in range(len(xs)):
                ax.plot(xs[i],ys[i],'o',label=labels[i],color=coloration[i])
        elif average == 'scatter':
            for i in range(len(xs)):
                ax.plot(xs[i],ys[i],'o',label=labels[i],color=coloration[i],markersize=2)
        else:
            for i in range(len(xs)):
                plt.errorbar(xs[i],ys[i],xerr=xerrs[i],yerr=yerrs[i],label=labels[i],\
                             color = coloration[i],fmt = fmt,capsize = 3)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()

    def plot_err(self,ion,ax=None,show=True,**kwargs):
        print("Current constraints (name): "+self.currentQuasarArrayName)
        ax = ax or plt.gca()
        plot_type = self.decide_plot_type(ion,**kwargs)
        ion,ion_name,labels = self.get_labels_from_ion(plot_type,ion,**kwargs)
        xarys,yarys = self.get_sightline_xy_vals(plot_type,ion,**kwargs)
        xs,ys,xerrs,yerrs,empty = self.process_datapoints(plot_type,ion,xarys,yarys,**kwargs)
        xlabel,ylabel,title = self.get_title_and_axislabels(plot_type,ion,ion_name,**kwargs)
        if not empty:
            self.plot_on_ax(ax,plot_type,xs,ys,labels,xerrs,yerrs,xlabel,ylabel,title,**kwargs)
            if show:
                plt.show()
        else:
            print "No values detected!"
    
    def sort_by_2D(self, criteria_x, criteria_y, bins_x = [0,np.inf], bins_y = [0,np.inf], \
                   reset = False, exploration_mode = False, atEnd_x = False, atEnd_y = False, \
                   onlyNonempty = False, splitEven_x = 0, splitEven_y = 0, reverse = False):
        
        #conduct safety checks on parameters
        bins_x, atEnd_x = self.prepare_plot(criteria_x, bins_x, atEnd_x)
        bins_y, atEnd_y = self.prepare_plot(criteria_y, bins_y, atEnd_y)
        
        
        #conduct a sort on one dimension and get the labels and bins for the other
        labels_y, bins_y, quasarbins_y = self.sort_by(criteria_y, bins_y, reset = reset,
                                                      exploration_mode = exploration_mode,
                                                      atEnd = atEnd_y, onlyNonempty = 
                                                      onlyNonempty,splitEven = splitEven_y,
                                                      reverse = False)
        labels_x, bins_x, _ = self.sort_by(criteria_x, bins_x, reset = reset,
                                                      exploration_mode = exploration_mode,
                                                      atEnd = atEnd_x, onlyNonempty = 
                                                      onlyNonempty,splitEven = splitEven_x,
                                                      reverse = False)
        quasarBins = np.zeros((len(bins_y)-1,len(bins_x)-1), dtype = object)  
                        

        #sort the other dimension and store sorted quasarbins
        for i in range(len(quasarbins_y)):
            quasarbin = quasarbins_y[i]      
            sorter = MultiSphereSorter(quasarbin, exploration_mode = exploration_mode)
            _,_, quasarbins_x = sorter.sort(criteria_x,bins_x,atEnd = atEnd_x)
           
            for j in range(len(quasarbins_x)):
                quasarBins[i][j] = quasarbins_x[j]
        
        #conduct post-sort checks
        quasarBins,bins_x,labels_x = self.postprocess_plot(quasarBins,bins_x,labels_x,onlyNonempty,reverse)
        quasarBins,bins_y,labels_y = self.postprocess_plot(quasarBins,bins_y,labels_y,onlyNonempty,reverse)
        
        return labels_x,labels_y,bins_x,bins_y,quasarBins
    
    
       
                              
    
    def faberplot(xVar,yVar,labels_x,labels_y,quasarArray):#...other args from plot_err or plot_hist):
        #after using sort_by_2d to get a set of labels and a 2d array of quasarspheres,
        #ask plot_err or plot_hist for completed plots of type given, for 
        #quasars in that cell of quasarArray, put them in subplots of a n by m subplots object
        #and show that plot
        return plot

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
                  save_fig = False, tolerance = 1e-5,ns = (42,15),logx = False, logy = True, plot_empties = False):
        hotcustom = self.definecolorbar()
        plt.figure()
        if isinstance(ion,tuple):
            ion_name = ion[1]
            ion = ion[0]
        else:
            ion_name = ion
        plt.register_cmap(cmap=hotcustom)
        if len(self.currentQuasarArray) == 0:
            print "no quasarspheres"
            return None
        if rlims is None:
            rlims = [0.1,np.inf]
        elif rlims == "all":
            rlims = [0.0,np.inf]
        vardict = {"theta":1,"phi":2,"r":3,"rdivR":3}
        if xVar == 'rdivR':
            self.constrain_current_Quasar_Array("Rvir_is_real",['True'],changeArrayName=False)
            if len(self.currentQuasarArray) == 0:
                print "You don't know Rvir for any galaxies! Plot with different xVar from 'rdivR'."
                return
        distances = "kpc" if xVar == "r" else "Rvir"
        gq = quasar_sphere.GeneralizedQuasarSphere(self.currentQuasarArray,distance = distances)
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
        xs = retxs
        self.debug.append(xs)
        
        if len(ys) == 0 and not plot_empties:
            print "No values found"
            return None
        H, xedges, yedges = np.histogram2d(xs, ys, bins=ns,weights = weight)
        H = H.T
        X, Y = np.meshgrid(xedges, yedges)
        plt.pcolormesh(X,Y, H, cmap=hotcustom)
        plt.colorbar(label = cbarlabel)
        ylabel,cd = self.get_ylabel_cd(ion,ion_name,islogy,1)
        plt.ylabel(ylabel)
        plt.xlabel(sightline_unit_labels[xVar])
        plt.title('%s Distribution %s'%(ion_name, extra_title))
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
    
    def prepare_plot(self, criteria, bins = [0,np.inf], reset = False, exploration_mode = False,\
        atEnd = False,onlyNonempty = False,splitEven = 0,reverse = False):
        if not (criteria in self.currentQuasarArray[0].__dict__.keys()):
            raise Exception("Criteria " + criteria + " does not exist. Please re-enter a valid criteria.")
    
        elif criteria == "ions":
            raise Exception("You cannot sort by 'ions'")
            
        if isinstance(bins, float) or isinstance(bins, int) or \
            (len(bins) == 1 and isinstance(bins[0], float)) or (len(bins) == 1 and isinstance(bins[0], int)):
            if isinstance(bins, list) or isinstance(bins,np.ndarray):
                bins = bins[0]
            bins = np.array([0.0, bins, np.inf])
        elif isinstance(bins, str) and criteria in stringcriteria:
            bins = [bins]
        if isinstance(atEnd,float):
            self.constrain_current_Quasar_Array("final_a0",[atEnd,np.inf],changeArrayName=False)
            if len(self.currentQuasarArray) == 0:
                print "No galaxies get to that high of a0"
            atEnd = True
        
        return bins, atEnd
    
    def postprocess_plot(self, quasarBins, bins = [0,np.inf], labels = None, reset = False, exploration_mode = False,\
        atEnd = False,onlyNonempty = False,splitEven = 0,reverse = False):
        if quasarBins is None:
            raise Exception("No quasars in quasarBin!")
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
            
        return labels, bins, quasarBins
    
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
                elif criteria != 'ions' and acceptedValues[i] == criteriaArray[j]:
                    add = True
                if add:
                    toAdd.append(self.array[j])
            resArray[i] = np.array(toAdd)
        return np.array(resArray)   
    
    def sort_by_default(self, criteria, bins, atEnd = False):
        if len(bins) == 1:
            criteriaArray = self.get_criteria_array(criteria, atEnd = atEnd)
            resArray = []
            for index in range(len(criteriaArray)):
                if criteriaArray[index] is np.nan:
                    continue
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
        criteriaArray = criteriaArray[criteriaArray>-1]
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
            min_a = min(final_a)
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
