import numpy as np
from quasarscan.plotting.read_and_init_from_file import read_and_init
from quasarscan.utils import utils
from quasarscan.utils.variable_lists import stringcriteria,intensives,intensiveslabels,\
                                            intensivespositions,sightline_xVars,param_xVars,\
                                            sightline_unit_labels,param_unit_labels,\
                                            all_known_variables
from quasarscan.plotting.sorter import MultiSphereSorter
from quasarscan.data_objects import gasbinning

class UnknownQTypeError(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)
class BadCriteriaError(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)
class SortFailureError(Exception):
    def __init__(self,message):
        self.message = message
        print(self.message)

class QuasarArrayHandler(object):
    def __init__(self,loadsim = "all",loadobs = 'all',loadempty = 'none'):
        self.sim_quasar_array = read_and_init(loadsim,'sim')
        self.current_sim_quasar_array = np.copy(self.sim_quasar_array)
        self.obs_quasar_array = read_and_init(loadobs,'obs')
        self.current_obs_quasar_array = np.copy(self.obs_quasar_array)
        self.empty_quasar_array = read_and_init(loadempty,'empty')
        self.current_empty_quasar_array = np.copy(self.empty_quasar_array)
        self.current_quasar_array_name = ''

    def get_current_quasar_array_name(self):
        return self.current_quasar_array_name if self.current_quasar_array_name!='' else 'all'

    def get_qlist(self,qtype):
        if qtype=='sim':
            return self.current_sim_quasar_array
        elif qtype=='obs':
            return self.current_obs_quasar_array
        elif qtype=='empty':
            return self.current_empty_quasar_array
        else:
            raise UnknownQTypeError('Tried to "get_qlist" for qtype: %s'%qtype)

    def get_qtype_index(self,qtype):
        try: 
            return {'sim':0,'obs':1,'empty':2}[qtype]
        except KeyError:
            raise UnknownQTypeError('Tried to "get_qtype_index" for qtype: %s'%qtype)

    def update_qlist(self,qtype,new_list,main_array=False):
        if not main_array:
            if qtype=='sim':
                self.current_sim_quasar_array = new_list
            elif qtype=='obs':
                self.current_obs_quasar_array = new_list
            elif qtype=='empty':
                self.current_empty_quasar_array = new_list
            else:
                raise UnknownQTypeError('Tried to "get_qtype_index" for qtype: %s'%qtype)
        else:
            if qtype=='sim':
                self.sim_quasar_array = new_list
            elif qtype=='obs':
                self.obs_quasar_array = new_list
            elif qtype=='empty':
                self.empty_quasar_array = new_list
            else:
                raise UnknownQTypeError('Tried to "get_qtype_index" for qtype: %s'%qtype)

    def length(self,include_nonsims):
        if include_nonsims:
            return (len(self.current_sim_quasar_array),len(self.current_obs_quasar_array),len(self.current_empty_quasar_array))
        return len(self.current_sim_quasar_array)

    #summary: print full list of QuasarSpheres, with metadata if requested (by 'criteria')
    #
    #inputs: *criteria: any criteria for a quasarSphere
    #
    #outputs: None, prints list of current quasarspheres
    def list_all_quasar_spheres(self, *criteria,qtype='sim',log=False):
        qlist = self.get_qlist(qtype)
        s = ""
        if len(qlist)==0:
            print("no QuasarSpheres")
        for q in sorted(qlist, key=lambda x: (x.fullname,x.rounded_redshift)):
            s += q.fullname +" at z=%s ("%q.rounded_redshift
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

    #summary: cancel all constraints
    #
    #inputs: None
    #        
    #outputs: None, changes state of mq.currentQuasarArray and mq.currentQuasarArrayName
    def reset_current_quasar_array(self):
        self.current_sim_quasar_array = np.copy(self.sim_quasar_array)
        self.current_obs_quasar_array = np.copy(self.obs_quasar_array)
        self.current_empty_quasar_array = np.copy(self.empty_quasar_array)
        self.current_quasar_array_name = ''

    #summary: Check if you can constrain by a parameter
    #
    #inputs: constrainCriteria:  what criteria to constrain by
    #        atEnd: list of lists of quasarSphere objects that fit in the bins
    #        
    #outputs: oldQuasarArray: copy of currentQuasarArray that can be modified without affecting original
    #         if constraint is illegal, raises error
    def check_criteria_works(self,constrain_criteria,qtype,at_end=False,**kwargs):
        qlist = self.get_qlist(qtype)
        if len(qlist)==0:
            #Cannot constrain further
            return
        elif constrain_criteria not in qlist[0].__dict__.keys():
            print("Constrain criteria " + constrain_criteria + " does not exist. Please re-enter a valid criteria.")
            raise BadCriteriaError(constrain_criteria)
        elif isinstance(at_end,float):
            old_qlist=np.copy(qlist)
            self.constrain_current_Quasar_Array("final_a0",bins=[at_end-.1,np.inf],change_array_name=False)
            if len(qlist) == 0:
                print("No galaxies get to that high of a0")
            self.update_qlist(qtype,old_qlist)
            return

    #summary: helper for 'sort_by'. Sets up bins if needed
    #
    #inputs: criteria: list of strings for labelling points in legend
    #        bins: the bins to compare to 
    #        atEnd: list of lists of quasarSphere objects that fit in the bins
    #        
    #outputs: bins: the bins to compare to 
    #         oldQuasarArray: copy of currentQuasarArray that can be modified without affecting original
    def prepare_to_sort(self, criteria, bins,qtype,at_end,**kwargs):
        
        self.check_criteria_works(criteria,qtype,at_end,**kwargs)
        if criteria == "ions":
            raise BadCriteriaError("You cannot sort by 'ions'")
            
        # can sort by giving one value, it will assume you mean "less than value" and "greater than value"
        if isinstance(bins, float) or isinstance(bins, int):
            bins = np.array([0.0, bins, np.inf])
        elif isinstance(bins, str) and criteria in stringcriteria:
            bins = [bins]
        return bins

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
    def postprocess_sorted(self, labels, bins, quasar_bins, reverse,only_nonempty = False,**kwargs):
        if quasar_bins is None:
            raise SortFailureError("No quasars in quasarBin!")
        empty = True
        nonempty_array = []
        nonempty_label_array = []
        for i,item in enumerate(quasar_bins):
            if len(item)>0:
                nonempty_array.append(item)
                nonempty_label_array.append(labels[i])
                empty = False
        if only_nonempty:
            labels, bins, quasar_bins = np.array(nonempty_label_array),bins, np.array(nonempty_array)
            print("All bins are empty." if empty else "")
        if reverse:
            labels = utils.reversearray(labels)
            bins = utils.reversearray(bins)
            quasar_bins = utils.reversearray(quasar_bins)
            
        return labels, bins, quasar_bins

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
    def sort_by(self, criteria, bins,at_end,split_even,reverse,sort_w_qtype = 'sim',**kwargs):
        ret_dict = {}

        default_list = ['sim','obs','empty']
        list_to_do = [default_list[default_list.index(sort_w_qtype)],\
                    default_list[(default_list.index(sort_w_qtype)+1)%3],\
                    default_list[(default_list.index(sort_w_qtype)+2)%3]]
        for i,qtype in enumerate(list_to_do):
            qlist = self.get_qlist(qtype)
            bins = self.prepare_to_sort(criteria,bins,qtype,at_end,**kwargs)
            sorter = MultiSphereSorter(qlist)
            if split_even and i==0:
                labels, bins, quasar_bins = sorter.split_even(criteria,split_even,at_end = at_end)
            elif split_even and i>0:
                _,_,quasar_bins = sorter.sort(criteria,bins,at_end = at_end)
            elif not split_even and i==0:
                labels, bins, quasar_bins = sorter.sort(criteria,bins,at_end = at_end)
            elif not split_even and i>0:
                _,_,quasar_bins = sorter.sort(criteria,bins,at_end = at_end)
            labels,bins,quasar_bins = self.postprocess_sorted(labels,bins,quasar_bins,reverse,**kwargs)
            ret_dict[qtype] = quasar_bins
        retlists = [ret_dict['sim'],ret_dict['obs'],ret_dict['empty']]
        return labels,bins,retlists,criteria

    def sort_by_2D(self, criteria_x,criteria_y, bins_x,bins_y,at_end_x,at_end_y,\
                    split_even_x,split_even_y,reverse_x,reverse_y,**kwargs):
        
        #conduct a sort on one dimension and get the labels and bins for the other
        labels_y, bins_y, quasarbins_y,_ = self.sort_by(criteria_y, 
                                                      bins_y,
                                                      at_end = at_end_y,
                                                      split_even = split_even_y,
                                                      reverse = ~reverse_y,**kwargs)
        labels_x, bins_x, _,_            = self.sort_by(criteria_x, 
                                                      bins_x, 
                                                      at_end = at_end_x,
                                                      split_even = split_even_x,
                                                      reverse = False)
        if reverse_x:
            labels_x = utils.reversearray(labels_x)
        quasarBins = [np.zeros((len(labels_y),len(labels_x)), dtype = object),
                      np.zeros((len(labels_y),len(labels_x)), dtype = object),
                      np.zeros((len(labels_y),len(labels_x)), dtype = object)]

        #sort the other dimension and store sorted quasarbins
        for i in range(3):
            for j,qlist_y in enumerate(quasarbins_y[i]):
                sorter = MultiSphereSorter(qlist_y)
                _,_, quasarbins_x = sorter.sort(criteria_x,bins_x,at_end = at_end_x)
                if reverse_x:
                    quasarbins_x = utils.reversearray(quasarbins_x)
                for k,qlist_x in enumerate(quasarbins_x):
                    quasarBins[i][j][k] = qlist_x
        return labels_x,labels_y,bins_x,bins_y,quasarBins  

    #summary: similar to 'constrain_current_Quasar_Array' but specific to gasbins
    #
    #inputs: gasbintype: only keep lines where this type of gasbin was defined
    #        
    #outputs: none, changes state of 'currentQuasarArray'    
    def constrain_via_gasbins(self,gasbintypes=None):
        if gasbintypes == []:
            return gasbintypes
        if gasbintypes == None:
            gasbintypes = input("Available bins are: %s. Please enter one."%gasbinning.possible_bin_types)
        if isinstance(gasbintypes,str):
            gasbintypes = [gasbintypes]
        if gasbintypes == []:
            return
        else:
            to_return = []
            for q in self.current_sim_quasar_array:
                for gasbintype in gasbintypes:
                    if gasbintype in q.gasbins:
                        to_return.append(q)
        if len(to_return)==0:
            print('all quasarspheres removed...')
        self.update_qlist('sim',np.array(to_return))
        
    #summary: helper function for constrain_array_helper
    #
    #inputs: constrainCriteria: what criteria to constrain by (simname, simnum, Mvir, Rvir, SFR, 
    #                           etc. See 'quasar_sphere.py' for full list)
    #        bins: either a list of two numbers, if a numerical criteria, or several strings if string criteria
    #              can leave as None if splitEven is used
    #        exclude: remove the ones listed, instead of keep the ones listed
    #        
    #outputs: bins: the bins to compare to 
    def get_bin_values(self,constrain_criteria,bins,qtype,exclude=False,**kwargs):
        if bins is None:
            bins = [0,np.inf]
        if isinstance(bins, list):                
            if len(bins) != 2 and constrain_criteria not in stringcriteria:
                raise BadCriteriaError("Length of bins must be 2: [lower,upper]")
        elif isinstance(bins,str) and constrain_criteria in stringcriteria:
            bins = [bins]
        if exclude:
            qlist = self.get_qlist(qtype)
            possible_bins = []
            for q in qlist:
                try:
                    possible_bins.append(eval("q.%s"%constrain_criteria))
                except AttributeError:
                    raise BadCriteriaError("qtype '%s' does not have property '%s'"%(qtype,constrain_criteria))
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
    def constrain_array_helper(self,sorters,constrain_criteria,bins,qtypes,split_even=False,at_end=False,\
                               set_main_array=False,**kwargs):
        for i,qtype in enumerate(qtypes):
            if split_even and i==0:
                labels,bins,temp = sorters[qtype].split_even(constrain_criteria,2,at_end=at_end)
                if split_even == "high":
                    take = 1
                elif split_even == "low":
                    take = 0
                else:
                    raise SortFailureError("please use split_even = 'low' or 'high'")
                labels = np.array([labels[take]])
                bins = np.array([bins[take],bins[take+1]])
                temp = np.array([temp[take]])
            else:
                labels, bins, temp = sorters[qtype].sort(constrain_criteria,bins,at_end=at_end)
            if len(temp)>0:
                update = np.concatenate(temp)
            else:
                update = np.array([])
            self.update_qlist(qtype,update)
            if set_main_array:
                self.update_qlist(qtype,update,main_array=True)
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
    def change_array_name(self,constrain_criteria,bins,change_array_name=True,exclude=False,**kwargs):
         if change_array_name:
            self.current_quasar_array_name += constrain_criteria 
            if not constrain_criteria in stringcriteria:
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
                    self.current_quasar_array_name += "any"
                else:
                    self.current_quasar_array_name += "%s%s"%(lowlabel,highlabel)
            else:
                if exclude:
                    self.current_quasar_array_name += 'exclude'
                for accepted_value in bins:
                    self.current_quasar_array_name += accepted_value.replace(" ","")

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
    def constrain_current_quasar_array(self, constrain_criteria,bins=None,qtype='all',**kwargs):
        if bins == []:
            return bins
        if isinstance(qtype,list):
            qtypes = qtype
        else:
            qtypes = ['sim','obs','empty'] if qtype == 'all' else [qtype]
        sorters = {}
        for qtype in qtypes:
            qlist = self.get_qlist(qtype)
            self.check_criteria_works(constrain_criteria,qtype,**kwargs)
            bins = self.get_bin_values(constrain_criteria,bins,qtype,**kwargs)
            sorters[qtype] = MultiSphereSorter(qlist)
        bins = self.constrain_array_helper(sorters,constrain_criteria,bins,qtypes,**kwargs)
        self.change_array_name(constrain_criteria,bins,**kwargs)
        return bins

    def impose_requirements(self,filter_for,qtype,verbose=False):
        if verbose and filter_for != []:
            print('Restricting quasars of type "%s" to ones containing %s'%(qtype,filter_for))
        ion_constraints = []
        gasbin_constraints = []
        metadata_constraints = []
        for item in filter_for:
            if ':' in item:
                ion_name = item.split(':')[0]
                gasbin_name = item.split(':')[1]
                if gasbin_name in ['cdens','fraction','eb']:
                    continue
                assert utils.string_represents_ion(ion_name)
                ion_constraints.append(ion_name)
                gasbin_constraints.append(gasbin_name)
            elif utils.string_represents_ion(item):
                ion_constraints.append(item)
            elif item in param_xVars:
                metadata_constraints.append(item)
        self.constrain_current_quasar_array('ions',ion_constraints,qtype=qtype)
        self.constrain_via_gasbins(gasbin_constraints)
        for item in metadata_constraints:
            self.constrain_current_quasar_array(item,qtype=qtype,change_array_name=False)

 




