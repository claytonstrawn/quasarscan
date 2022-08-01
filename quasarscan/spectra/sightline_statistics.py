#for one code:
#iterate over all the lines
#if any lines have the ion
#count how many components within that line
#but only counting ones where this ion is found
#return the average number

def comps_fraction(ion,comp_list):
    line_comp_tot = 0
    comp_fraction = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                line_comp_tot = line_comp_tot + 1
                break
    if len(comp_list) != 0:
        comp_fraction = line_comp_tot/len(comp_list)
    return comp_fraction

def ncomps(ion,comp_list):
    line_comp_tot = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                line_comp_tot = line_comp_tot + 1
                break
    return line_comp_tot

#if any lines have the ion
#collect b for each component within that line
#but only counting ones where this ion is found
#return the average b
def bparam(ion, comp_list):
    total_b = 0
    line_b_avg = 0
    counter = 0
    for comp in comp_list:
        for line in comp.list_of_lines:
            if ion == line.ion:
                counter += 1
                total_b = total_b + line.b
    if counter != 0:
        line_b_avg = total_b / counter
    return line_b_avg