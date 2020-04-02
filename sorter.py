import numpy as np

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


class MultiSphereSorter(object):
    def __init__(self,myArray):
        self.array = myArray
    #param bins the lower parameter is inclusive, the upper parameter is exclusive
    def sort(self, criteria, bins,atEnd = False):
        labels = self.make_labels(criteria, bins,atEnd = atEnd)
        if criteria in stringcriteria:
            sortfn = self.sort_by_strparam
        else:
            sortfn = self.sort_by_default
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
            booleanindices = np.logical_and(criteriaArray >= float(bins[index]), criteriaArray < bins[index+1]) 
            toadd = self.array[booleanindices]
            resultBins[index] = toadd
        return np.array(resultBins)
    
    def splitEven(self,criteria,num,atEnd = False):
        if criteria in stringcriteria:
            print("cannot splitEven over string criteria")
        criteriaArray = self.get_criteria_array(criteria, atEnd = atEnd)
        # criteriaArray = criteriaArray[criteriaArray>-1]
        sortedcriteriaArray = np.sort(criteriaArray)
        if atEnd:
            sortedcriteriaArray = np.unique(sortedcriteriaArray)
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
