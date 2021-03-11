import numpy as np

from quasarscan.utils.variable_lists import stringcriteria,intensives,intensiveslabels,\
                                            intensivespositions,sightline_xVars,param_xVars,\
                                            sightline_unit_labels,param_unit_labels,\
                                            all_known_variables

class IllegalSort(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)

class MultiSphereSorter(object):
    def __init__(self,my_array):
        self.array = my_array
    #param bins the lower parameter is inclusive, the upper parameter is exclusive
    def sort(self, criteria, bins,at_end = False):

        labels = self.make_labels(criteria, bins,at_end = at_end)
        if criteria in stringcriteria:
            sortfn = self.sort_by_strparam
        else:
            sortfn = self.sort_by_numparam
        return labels, bins, sortfn(criteria, bins, at_end = at_end)

    def sort_by_strparam(self, criteria, accepted_values,at_end = False):
        if at_end:
            raise IllegalSort('Cannot use at_end with string-based criteria %s'%criteria)
        criteria_array = self.get_criteria_array(criteria)
        result_array = np.empty(len(accepted_values),dtype = 'object')
        for i,value in enumerate(accepted_values):
            to_add = []
            for j,qs_values in enumerate(criteria_array):
                if criteria == "ions" and value in qs_values:
                    to_add.append(self.array[j])
                elif criteria != 'ions' and value == qs_values:
                    to_add.append(self.array[j])
            result_array[i] = np.array(to_add)
        return np.array(result_array,dtype=object)   
    
    def sort_by_numparam(self, criteria, bins, at_end = False):
        result_bins = [None] * (len(bins)-1)
        criteria_array = self.get_criteria_array(criteria,at_end=at_end)
        for index in range(len(bins)-1):
            booleanindices = np.logical_and(criteria_array >= float(bins[index]), criteria_array < bins[index+1]) 
            to_add = self.array[booleanindices]
            result_bins[index] = to_add
        return np.array(result_bins,dtype=object)
    
    def split_even(self,criteria,num,at_end = False,verbose = False):
        if criteria in stringcriteria:
            print("cannot splitEven over string criteria")
        criteria_array = self.get_criteria_array(criteria, at_end = at_end)
        # criteriaArray = criteriaArray[criteriaArray>-1]
        sorted_criteria_array = np.sort(criteria_array)
        if at_end:
            sorted_criteria_array = np.unique(sorted_criteria_array)
        quotient = len(sorted_criteria_array) // num
        if quotient == 0 and verbose:
            print("Warning: Number of bins exceeds length of criteria array. Not all bins will be filled.")
        remainder = len(sorted_criteria_array) % num
        bin_edges = np.zeros(num + 1)
        bin_edges[0] = 0
        bin_edges[-1] = np.inf
        j = 0
        for i in range(1,num):
            if i <= remainder:
                bin_edges[i] = np.mean(sorted_criteria_array[i*quotient+j : i*quotient+2+j])
                j+=1
            else:
                bin_edges[i] = np.mean(sorted_criteria_array[i*quotient-1+j : i*quotient+1+j])
        return self.sort(criteria, bin_edges,at_end = at_end)
    
    
    def get_criteria_array(self, criteria,at_end = False):
        if not at_end:
            res = []
            for q in self.array:
                criteria_vals = eval("q." + criteria)
                res.append(criteria_vals)
            return np.array(res,dtype=object)
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
    
    def make_labels(self,criteria, bins, at_end=False):
        labels = []
        if criteria in stringcriteria:
            labels = bins
        else:
            if at_end:
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
                    unique_name = "%s < %s"%(criteria, highstr)
                elif high == np.inf:
                    unique_name = "%s > %s"%(criteria, lowstr)
                else:
                    unique_name = "%s < %s < %s"%(lowstr,criteria,highstr)
                labels.append(unique_name)
        return labels
