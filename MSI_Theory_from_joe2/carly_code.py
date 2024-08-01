#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import cantera as ct

class MSI_Theory:
    
    def __init__(self, calculation_directory, p_list, states_dict, TP_shape, TP_dict):
        self.calculation_directory = calculation_directory
        self.p_list = p_list
        self.states_dict = states_dict  
        self.TP_shape = TP_shape
        self.TP_dict = TP_dict    
    
    def extract_rate_constants(self):
        
        def load_text_file(self, file_name):
            with open(file_name) as f:
                lines = [line.rstrip() for line in f]
            return lines

        def pull_out_chebyshev_fit(self, lines,p_list,shape):
            chebyshev_fit=[]
            reaction_label=[]
            counter=None
            for i,l in enumerate(lines):
                if l.strip() in p_list:
                    reaction_label.append(l.strip())
                    temp = lines[i+2]
                    temp= temp.replace('[[','')
                    temp= temp.replace(']]','')
                    l_split = temp.split(' ')
                    temp_list = []
                    for i,nmb in enumerate(l_split):
                        if nmb != '':
                            nmb = nmb.replace(',','')
                            temp_list.append(float(nmb))
                    
                    arr = np.array(temp_list)
                    arr = arr.reshape(shape)
                    chebyshev_fit.append(arr)

            return reaction_label,chebyshev_fit

        def build_reaction_string(self, p_list,temp_dict):
            reactions=[]
            for string in p_list:
                split = string.split('->')
                reactant = temp_dict[split[0]]
                product = temp_dict[split[1]]
                final = reactant +' <=> ' + product
                reactions.append(final)
            return reactions

        def convert_units(self, dict_of_reactions,constant = np.log10(6.0221409e+23)):
            for key in dict_of_reactions.keys():
                temp = key.split('<=>')
                #print(temp[0][0])
                if '+' in temp[0] or temp[0][0]=='2':
                    print(key)
                    cheby_fit = dict_of_reactions[key]
                    cheby_fit[0,0] = cheby_fit[0,0]+constant
                    dict_of_reactions[key] = cheby_fit   
                
            return dict_of_reactions

        def write_out_reactions(self, dict_of_reactions_unites_converted,T_min,T_max,P_min,P_max):
            with open('readme.txt', 'w') as f:
                for i, reaction in enumerate(dict_of_reactions_unites_converted.keys()):
                    f.write('#  Reaction '+str(i)+'\n')
                    f.write('chebyshev_reaction(\''+reaction+'\',\n')
                    f.write('T_min='+'%s' % float('%.5g' % T_min)+', T_max='+'%s' % float('%.5g' % T_max)+',\n')
                    f.write('P_min=('+'%s' % float('%.5g' % P_min)+', \'atm\'),')
                    f.write('P_max=('+'%s' % float('%.5g' % P_max)+', \'atm\'),\n')
                    f.write('coeffs=[')
                    arry_temp = dict_of_reactions_unites_converted[reaction]
                    tempvar=np.arange(len(arry_temp[:,0]))
                    for i in np.arange(len(arry_temp[:,0])):
                        if i==0:
                            converted_first_row = arry_temp[i,:]
                            if converted_first_row.ndim ==1:
                                converted_first_row[0] = converted_first_row[0]+0
                            else:
                                converted_first_row[0,0] = converted_first_row[0,0] +0
                            f.write(str(list(converted_first_row)))
                        elif i!=0 and i!=tempvar[-1]:
                                f.write('                           '+str(list(arry_temp[i,:])))
                        if i!=tempvar[-1]:
                                f.write(',\n')
                        elif i==tempvar[-1]:
                                f.write('                           '+str(list(arry_temp[i,:])))
                                f.write('])\n\n')
                    f.write('\n')
                    f.write('\n')
                    
        # def run(self, calculation_directory, p_list, states_dict, TP_shape, TP_dict):
        file_name = 'Chebyshev_fit.txt'
        lines=load_text_file(self.calculation_directory + file_name)
        reaction_label, fit = self.pull_out_chebyshev_fit(lines, self.p_list, self.TP_shape)
        reactions = self.build_reaction_string(self.p_list, self.states_dict)
        dict_of_reactions = dict(zip(reactions,fit))
        dict_of_reactions_unites_converted = self.convert_units(dict_of_reactions)

        write_out_reactions(dict_of_reactions_unites_converted, T_min=200, T_max=4000, P_min=0.0001, P_max=100)

    def sensitivity_parse(self):

        def load_text_file(self, file_name):
            with open(file_name) as f:
                lines = [line.rstrip() for line in f]
            return lines

        def build_dictonary_by_key_word_paramter(self, key_words,number_of_sets,lines):
            key_word_stored_sens_lists=[]
            counter=None
            for l in lines:
                if l in key_words:
                    counter = 0
                    key_word_temp_list = []
                if '[[' in l and ']]' in l:
                    counter+=1
                    temp_list = []
                    appnd = True
                    l = l.replace('[[','')
                    l = l.replace(']]','')
                    l_split = l.split(' ')
                    for i,nmb in enumerate(l_split):
                        if nmb != '':
                            nmb = nmb.replace(',','')
                            temp_list.append(float(nmb))
                    key_word_temp_list.append(temp_list)
                    #counter+=1
                    
                if counter == number_of_sets:
                    key_word_stored_sens_lists.append(key_word_temp_list)
                    counter=0
            
            dict_by_key_words = dict(zip(key_words,key_word_stored_sens_lists))
            return dict_by_key_words

        def build_temp_dict(self, key_words):
            temp_dict = {}
            for word in key_words:
                temp_dict[word] = []
            
            return temp_dict

        def convert_units(self, dict_by_key_words,temp_dict):
            for key in dict_by_key_words:
                for i,reaction_set in enumerate(dict_by_key_words[key]):
                    if 'energy' in key or 'Energy' in key:
                        temp_arr = np.array(reaction_set)
                        temp_arr = temp_arr.reshape((temp_arr.shape[0],1))
                        zeros = np.zeros((temp_arr.shape[0],1))
                        temp_arr = np.hstack((temp_arr,zeros))
                        temp_arr = temp_arr* np.log(10)* 349.757
                        temp_dict[key].append(temp_arr)
                    else:
                        temp_arr = np.array(reaction_set)
                        temp_arr = temp_arr.reshape((temp_arr.shape[0],1))
                        zeros = np.zeros((temp_arr.shape[0],1))
                        temp_arr = np.hstack((temp_arr,zeros))
                        temp_arr = temp_arr* np.log(10)
                        temp_dict[key].append(temp_arr)

            return temp_dict

        def convert_units_two(self, dict_by_key_words,temp_dict,shape=(6,3)):
            for key in dict_by_key_words:
                for i,reaction_set in enumerate(dict_by_key_words[key]):
                    if 'energy' in key or 'Energy' in key:
                        temp_arr = np.array(reaction_set)
                        temp_arr = temp_arr.reshape(shape)
                        temp_arr = temp_arr* np.log(10)* 349.757
                        temp_dict[key].append(temp_arr)
                    else:
                        temp_arr = np.array(reaction_set)
                        temp_arr = temp_arr.reshape(shape)
                        temp_arr = temp_arr* np.log(10)
                        temp_dict[key].append(temp_arr)

            return temp_dict

        def build_empty_nested_list(self, temp_dict,key_words):
            empty_nested_list =[]
            for nmb_reactions in range(len(temp_dict[key_words[0]])):
                empty_nested_list.append([])
            return empty_nested_list
            
            
        def breakup_by_reaction(self, empty_nested_list,key_words,temp_dict):
            for i in range(len(empty_nested_list)):
                for word in key_words:
                    empty_nested_list[i].append(temp_dict[word][i])
            return empty_nested_list

        def pull_out_key_words(self, lines):
            key_word=[]
            counter=None
            for i,l in enumerate(lines):
                if '====' in l:
                    key_word.append(lines[i+1])

            return key_word

        def build_reaction_string(self, p_list,temp_dict):
            reactions=[]
            
            for string in p_list:
                split = string.split('->')
                print(split)
                reactant = temp_dict[split[0]]
                product = temp_dict[split[1]]
                final = reactant +' <=> ' + product
                reactions.append(final)
            return reactions

        def reduced_T(self, T, T_min, T_max):
            '''Calculate the reduced temperature.'''
            T = np.array(T)
            T_tilde = 2. * T ** (-1.0) - T_min ** (-1.0) - T_max ** (-1.0)
            T_tilde /= (T_max ** (-1.0) - T_min ** (-1.0))
            return T_tilde

        def calc_reduced_P(self, P,P_min,P_max):
                
            numerator = 2*np.log10(P) - np.log10(P_min) - np.log10(P_max)
            denominator = np.log10(P_max) - np.log10(P_min)
            P_reduced = np.divide(numerator,denominator)
            return P_reduced

        def calc_chevy(self, T,P,alpha):
            #calculate rate constants helper function
            T_reduced_list = self.reduced_T(T,200,2400)
            P_reduced_list = self.calc_reduced_P(P,.0001,10.0)
            values = np.polynomial.chebyshev.chebval2d(T_reduced_list,P_reduced_list,alpha)

            return values

        def plot_sens(self, dictonary,path,key_words):
            for key in dictonary.keys():
                for j, sens in enumerate(dictonary[key]):
                    plt.figure()
                    values = self.calc_chevy(np.arange(500,2400),[1]*(2400-500),sens)
                    temp_list = np.arange(500,2400)
                    #p_list = np.arange(.0001,10.0)
                    plt.plot(temp_list,values)
                    plt.title('Reaction: '+key+' '+ 'Parameter: '+ key_words[j])
                    plt.xlabel('Temperature (K)')
                    plt.ylabel('dln(k)/dln(paramter)')
                    plt.savefig(path+'/'+key+'_'+key_words[j]+'.pdf',bbox_inches='tight')
                    
        def plot_sens2(self, dictonary,key_words,pressure_list, save_fig = False):
            for key in dictonary.keys():
                for j, sens in enumerate(dictonary[key]):
                    plt.figure()
                    for pressure in pressure_list:          
                        values = self.calc_chevy(np.arange(290,3000),[pressure]*2710,sens)
                        temp_list = np.arange(290,3000)
                        plt.plot(temp_list,values,label=str(pressure)+' atm')
                        plt.legend()
                        
                    plt.title('Reaction: '+key+' '+ 'Parameter: '+ key_words[j])
                    plt.xlabel('Temperature (K)')
                    plt.ylabel('dln(k)/dln(paramter)')
                    if save_fig == True:
                        plt.savefig(key+'_'+key_words[j]+'.pdf',bbox_inches='tight')
                    
        file_name = 'Chebyshev_sens.txt'
        lines=self.load_text_file(self.calculation_directory+file_name)
        self.key_words = self.pull_out_key_words(lines)
        reaction_list = self.build_reaction_string(self.p_list, self.states_dict) #{'P1':'N2O + O', 'P2':'N2 + O2'}
        number_of_sets = len(self.p_list)
        dict_by_key_words = self.build_dictonary_by_key_word_paramter(self.key_words,number_of_sets,lines)
        temp_dict = self.build_temp_dict(self.key_words) 
        temp_dict = self.convert_units(dict_by_key_words,temp_dict) # 
        # temp_dict = convert_units_two(dict_by_key_words,temp_dict,shape=(7,1))

        empty_nested_list = self.build_empty_nested_list(temp_dict, self.key_words)

        breakup_by_reaction = self.breakup_by_reaction(empty_nested_list, self.key_words, temp_dict)

        self.test = dict(zip(reaction_list, breakup_by_reaction))

        return self.key_words, self.test

    def vibe_checker(self, reaction, channel, original_model, chebyshev_model, temperature_list = np.arange(200,4000), pressure_list = [1], conversion = 1000):

        plot_chevy_fit = False
        plot_arrhenius_fit = True
        plot_mess_data = True

        gas_chebyshev = ct.Solution(chebyshev_model)
        gas_original = ct.Solution(original_model)

        gas_chebyshev.reaction_equations()

        cheby_index = gas_chebyshev.reaction_equations().index(reaction)
        k=[]
        for i,p in enumerate(pressure_list):
            k_temp = []
            for j,temperature in enumerate(temperature_list):
                gas_chebyshev.TPX = temperature,p*101325,{'Ar':1}
                k_temp.append(gas_chebyshev.forward_rate_constants[cheby_index]*conversion)
            k.append(k_temp)
            
        regular_index = gas_original.reaction_equations().index(reaction)
        k2=[]
        for i,p in enumerate(pressure_list):
            k_temp2 = []
            for j,temperature in enumerate(temperature_list):
                gas_original.TPX = temperature,p*101325,{'Ar':1}
                k_temp2.append(gas_original.forward_rate_constants[regular_index]*conversion)
            k2.append(k_temp2)    
            
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=(7,7))
        fig.subplots_adjust(hspace=0)

        color_list = ['b','k','r','g','pink']
                
        for i,p in enumerate(pressure_list):
            
            axs[0].semilogy(temperature_list,k2[i],color='r',label='original')
            axs[0].semilogy(temperature_list,k[i],'--',color='b',label='chebyshev')
            # plt.xlim(1000,2500)
            # plt.ylim(.01,10000000)
            axs[0].set_title(reaction)
            axs[0].set_xlabel('Temperature [K]')
            axs[0].set_ylabel('k')
            axs[0].tick_params(axis = 'x', direction='in')
            axs[0].tick_params(axis = 'y', direction='in')
            

            
        df = pd.read_csv(os.path.join(self.calcuation_directory,'T_rate.csv'), skiprows=3)
        # df2 = df.tail(df.shape[0] -3)

        df_temp = df['T']
        df_k = df[channel]

        df_k_converted = np.multiply(df_k,6.022e23)

        axs[0].scatter(df_temp, df_k_converted,label='raw data',color='k')
        axs[0].legend()


        # plt.figure()    
        for i,p in enumerate(pressure_list):
            
            
            k_reshaped = np.interp(df_temp, temperature_list, k[i])
            percent_change = 100 * np.divide(np.subtract(df_k_converted, k_reshaped), df_k_converted)  
            # difference = np.subtract(k_reshaped,df_k_converted)
            
            k2_reshaped = np.interp(df_temp, temperature_list, k2[i])
            percent_change2 = 100 * np.divide(np.subtract(df_k_converted, k2_reshaped), df_k_converted)         
            
            
            # axs[0].title(reaction + ' Percent Change')
            axs[1].plot(df_temp,percent_change, 'b--')
            axs[1].plot(df_temp,percent_change2, 'r-')
            axs[1].set_xlabel('Temperature [K]')
            axs[1].set_ylabel('Percent Change')
            axs[1].tick_params(axis = 'x', direction='in')
            axs[1].tick_params(axis = 'y', direction='in')
            
        fig.savefig('vibe_check.pdf',dpi=1000,bbox_inches='tight')
        
        
    def pressure_checker(self, reaction, channel, original_model, chebyshev_model, temperature_list = np.arange(200,4000), pressure_list = [1], conversion = 1000):

        plot_chevy_fit = False
        plot_arrhenius_fit = True
        plot_mess_data = True

        gas_chebyshev = ct.Solution(chebyshev_model)
        gas_original = ct.Solution(original_model)

        gas_chebyshev.reaction_equations()

        reaction = 'N2O + O <=> 2 NO'

        temperature_list = np.arange(200,4000)
        pressure_list = [1,5,20]
        lines = ['solid', (0, (5, 10)), 'dotted']
        colors = ['r','b','g','k']
        cheby_index = gas_chebyshev.reaction_equations().index(reaction)
        k=[]
        for i,p in enumerate(pressure_list):
            k_temp = []
            for j,temperature in enumerate(temperature_list):
                gas_chebyshev.TPX = temperature,p*101325,{'Ar':1}
                k_temp.append(gas_chebyshev.forward_rate_constants[cheby_index]*1000)
            k.append(k_temp)
            
        regular_index = gas_original.reaction_equations().index(reaction)
        k2=[]
        for i,p in enumerate(pressure_list):
            k_temp2 = []
            for j,temperature in enumerate(temperature_list):
                gas_original.TPX = temperature,p*101325,{'Ar':1}
                k_temp2.append(gas_original.forward_rate_constants[regular_index]*1000)
            k2.append(k_temp2)    
            
        fig, axs = plt.subplots(2, 1, sharex=True, figsize=(7,7))
        fig.subplots_adjust(hspace=0)


                
        # for i,p in enumerate(pressure_list):
            
        #     axs[0].semilogy(temperature_list,k2[i],linestyle = lines[i],color='r',label='gonzalez ' + str(p) + ' atm')
        #     axs[0].semilogy(temperature_list,k[i],linestyle = lines[i],color='b',label='mess ' + str(p) + ' atm')
            # plt.xlim(1000,2500)
            # plt.ylim(.01,10000000)
        axs[0].set_title(reaction)
        axs[0].set_xlabel('Temperature [K]')
        axs[0].set_ylabel('k [cm^3/mol*s]')
        axs[0].tick_params(axis = 'x', direction='in')
        axs[0].tick_params(axis = 'y', direction='in')
            

            
        df = pd.read_csv('/home/jl/MSI_Theory/PAPR-MESS_calculation/NO+NO/calculation_9/nominal/T_rate.csv', skiprows=3)

        k_lists = []
        T_list = []
        k_list = []
        signal = -4
        for i, val in enumerate(df.iloc[:,0]):
            
            if val == "========================================":
                signal = i
                k_lists.append([T_list,k_list])
                T_list = []
                k_list = []
                pass
            elif i - signal <= 3:
                pass
            else:
                T_list.append(eval(val))
                k_list.append(eval(df.iloc[i,1]))
        k_lists.append([T_list,k_list])
                
        # df_temp = df['T']
        # df_k = df['P1->P2']

        for i,p in enumerate(pressure_list):

            df_k_converted = np.multiply(k_lists[i][1],6.022e23)

            df_k_converted_0 = np.multiply(k_lists[0][1],6.022e23)

                
                
            # axs[0].semilogy(temperature_list,k[i],linestyle = lines[i],color='b',label='mess ' + str(p) + ' atm')
            

            axs[0].semilogy(k_lists[i][0], df_k_converted, linestyle = lines[i], label='raw data ' + str(p) + ' atm', color=colors[i])
            


            # plt.figure()    
            # for i,p in enumerate(pressure_list):
                
            if i != 0:
                
                # k_reshaped = np.interp(k_lists[i][0], temperature_list, k[i])
                percent_change = 100 * np.divide(np.subtract(df_k_converted, df_k_converted_0), df_k_converted)  
            # difference = np.subtract(k_reshaped,df_k_converted)
            
            # k2_reshaped = np.interp(df_temp, temperature_list, k2[i])
            # percent_change2 = 100 * np.divide(np.subtract(df_k_converted, k2_reshaped), df_k_converted)         
            
            
            # axs[0].title(reaction + ' Percent Change')
            # axs[1].plot(df_temp,percent_change, 'b--')
                axs[1].plot(k_lists[i][0],percent_change, linestyle = lines[i],  color=colors[i])
                
                

        axs[0].semilogy(temperature_list,k2[0],linestyle = 'dashed',color='k',label='gonzalez')
        k2_reshaped = np.interp(k_lists[0][0], temperature_list, k2[i])
        percent_change2 = 100 * np.divide(np.subtract(k2_reshaped, df_k_converted_0), k2_reshaped)  
        axs[1].plot(k_lists[0][0],percent_change2, linestyle = 'dashed',  color='k')
                

        axs[0].legend()
        axs[1].set_xlabel('Temperature [K]')
        axs[1].set_ylabel('Percent Change')
        axs[1].tick_params(axis = 'x', direction='in')
        axs[1].tick_params(axis = 'y', direction='in')
            
        # fig.savefig('Pressure_Dependence_Checker_N2O+O_N2+O2.pdf',dpi=1000,bbox_inches='tight')