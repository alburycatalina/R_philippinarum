# Original paper: Flye-Sainte-Marie, J., Jean, F., Paillard, C., Ford, S., Powell, E., Hofmann, E., Klinck, J., 2007. Ecophysiological dynamic model of individual growth of Ruditapes philippinarum. Aquaculture 266, 130–143. https://doi.org/10.1016/j.aquaculture.2007.02.017

# Script contributors: 
    # Catalina Albury (MSc Candidate, Dept. of Biology, Dalhousie University)
    # Diego Ibarra (Diego Ibarra Post-Doctoral Fellow, Dept. of Oceanography, Dalhousie University)

#%%
def load_defaults():

    print('Loading defaults...')

    # Framework

    days = 365 * 7 # Two years

    dt   = 0.001 # units: days



    # Parameters

    par = {}

    par['Filt_temp_16'] = 1.16 # max filtration at 16° C

    par['a_f'] = 20.049 # Allometric coefficient

    par['b_f'] = 0.257 # Allometric exponent

    par['AE_0'] = 0.1 # minimum additive assimilation rate

    par['AE_1'] = 0.6 # Max addtitive assimilation rate

    par['K_A'] = 2.74 # Half sat

    par['Resp_T_21'] = 0.24 # Maximum RT functiion @ 21° C

    par['a_r'] = 27.88 # Allometric coefficient of the equation relating respiration rate to DW

    par['b_r'] = 0.85 # Allometric exponent of the equation relating respiration rate to DW

    par['a_o'] = 851.88E-9 # Allometric coefficient of the equation relating average DW (W_o) to length

    par['b_o'] = 3.728 # Allometric coefficient of the equation relating average DW (W_o) to length

    par['a_om'] = 1703.76E-9 # Allometric coefficient of the equation relating maximum DW (W_max) to length

    par['L_max'] = 60 # Max clam length

    par['dl_o'] = 0.06 #0.0045 # Maximum length growth rate
        # This unit was analytically forced to make curve more realistic

    par['K_l'] = 0.09 #0.123 # Half-saturation constant of the Michaelis–Menten equation relating C to length growth rate

    par['C_s'] = -0.45 # Minimum value of C

    par['LengthRepMin'] = 20 # Minimum clam length for reproduction

    par['TempRepMin'] = 12 # Minimum temperature for reproduction

    par['ParGSI'] = 0.30 # GSI threshold value for partial spawning

    par['ParIC'] = 0.50 # Condition index threshold value for partial spawning

    par['ParPostIC'] = 0.23 # Condition index value after a partial spawning event

    par['SpawnRatio'] = 0.42 # GSI threshold value for principal spawning

    par['C2F'] = 0.0003 #0.2552 # Chla to dry ingestible organic matter conversion coefficient at Nole Station
        # We believe units were an issue for this parameter. Chlorophyll conversion factors that made sense where 10E-3 smaller than those presented

    par['RepEff_max'] = 0.5 # Guessed at what this could be. Not included in the paper as far as I can tell?
    


    # Initial conditions

    InitCond = {}

    InitCond['Weight'] = 0.1 # g DW

    InitCond['Length'] = 28 # mm

    InitCond['Condition_index'] = 50 # Unitless

    InitCond['Gonad'] = 0 # g DW

    InitCond['Temp'] = 15.0

    InitCond['chla'] = 1#12

    InitCond['RepEff'] = 0

    InitCond['GSI'] = 0
    
    return  days, dt, par, InitCond

#%%
def run(days, dt, par, InitCond):

    print('Running model...')

    # Import libraries

    import numpy as np

    import math



    # Setup the framework

    NoSTEPS = int(days / dt) # Calculates the number of steps

    time = np.linspace(0,days,NoSTEPS) # Makes vector array of equally spaced numbers



    # Create arrays of zeros

    Temp = np.zeros((NoSTEPS), float)

    Weight = np.zeros((NoSTEPS), float)

    Length = np.zeros((NoSTEPS), float)

    Condition_index = np.zeros((NoSTEPS), float)

    Gonad = np.zeros((NoSTEPS), float)

    RepEff = np.zeros((NoSTEPS), float)

    GSI = np.zeros((NoSTEPS), float)




    # Initializing with initial conditions

    Weight[0] = InitCond['Weight']

    Length[0] = InitCond['Length']

    Condition_index[0] = InitCond['Condition_index']

    Gonad[0] = InitCond['Gonad']

    RepEff[0] = InitCond['RepEff']

    GSI[0] = InitCond['GSI']

    Temp[0] = InitCond['Temp']



    for i in range(1, days):

        Temp[i] = (((20-10)/2) + 10) + ((20-10)/2)*np.sin((2*math.pi*i)/365)

        Temp[i] = 15





    # *****************************************************************************

    # MAIN MODEL LOOP *************************************************************

    for t in range(0,NoSTEPS-1):

        # Estimate Time rate of change of all State Variables

        # Pnet
        Filt_temp = -5.62e-3 * InitCond['Temp']**2 + 0.18 * InitCond['Temp'] - 0.30 # EQ 2
        
        Filt_weight = par['a_f'] * (Weight[t]**par['b_f']) # EQ 3
        
        Filt = (Filt_weight * (Filt_temp/par['Filt_temp_16'])) * 1.44  # EQ 4 (original)

        Food = InitCond['chla'] * par['C2F']

        Ingest = Food * Filt # EQ 5

        # AE = par['AE_0'] + ( (par['AE_1']* Weight[t]) / (par['K_A'] + Weight[t]) ) # EQ 6

        AE = 0.4 # EQ 6 (Catalina) - choose an AE between max and minimum because equation above is not working

        Assim = AE * Ingest # EQ 7
        
        Resp_temp = -4.75E-4 * Temp[t]**2 + 2.02E-2 * Temp[t] + 2.43E-2 # EQ 8

        Resp_weight = par['a_r'] * Weight[t]**par['b_r'] # EQ 9

        Resp = ((Resp_weight * Resp_temp)/par['Resp_T_21']) * 2.177E-3 # EQ 10
        
        # Removed duplicate inclusion of Filtration
        Pnet =  Assim - Resp # Eq 1 (Diego)
        
        
        # Condition Index :-)

        W_o = par['a_o'] * Length[t]**par['b_o'] # Eq 12

        W_max = par['a_om'] * Length[t]**par['b_o'] # Eq 13

        C = (Weight[t] - W_o) / (W_max - W_o) # Eq 11

        dCondition_indexdt = (C + 1.0562) / 0.0124


        # Length

        dLengthdt = par['dl_o'] *( (par['L_max'] - Length[t]) / par['L_max']) * ((C- par['C_s']) / (par['K_l'] + C - par['C_s'])) # Equation 14

        # Update and step (time-stepping) ------------------------------

        Weight[t+1] = Weight[t] + (Pnet * dt)

        Length[t+1] = Length[t] + (dLengthdt * dt)
        
        Condition_index[t+1] = Condition_index[t] + (dCondition_indexdt * dt)


        # Conservation of Mass Test
        # To test conservation of mass, subtract inputs (food) from outputs (respiration)
        # This curve should be equal to pnet
        
        conmass = Weight


    # end of main model LOOP***************************************************

    # *************************************************************************



    # Pack output into dictionary

    output = {}

    output['time'] = time

    output['Weight'] = Weight

    output['Length'] = Length

    output['Condition_index'] = Condition_index
    
    output['conmass'] = conmass



    print('Model run: DONE!!!')

    return  output

#%%





def plot_weight(output):

    import matplotlib.pyplot as plt

    from matplotlib import style

    # Plotting

    style.use('ggplot')

    fig, (ax) = plt.subplots(1,1)

    ax.plot(output['time']/365,output['Weight'],'r-')

    ax.set_ylabel('Weight (g DW)')

    ax.set_xlabel('Time (years)')

    plt.show()

    return



def plot_length(output):

    import matplotlib.pyplot as plt

    from matplotlib import style

    # Plotting

    style.use('ggplot')

    fig, (ax) = plt.subplots(1,1)

    ax.plot(output['time']/365,output['Length'],'b-')

    ax.set_ylabel('Length (mm)')

    ax.set_xlabel('Time (years)')

    plt.show()

    return



def plot_condition(output):

    import matplotlib.pyplot as plt

    from matplotlib import style

    # Plotting

    style.use('ggplot')

    fig, (ax) = plt.subplots(1,1)

    ax.plot(output['time']/365,output['Condition_index'],'g-')

    ax.set_ylabel('Condition Index (Unitless)')

    ax.set_xlabel('Time (years)')

    plt.show()

    return

def plot_masscon(output):

    import matplotlib.pyplot as plt

    from matplotlib import style

    # Plotting

    style.use('ggplot')

    fig, (ax) = plt.subplots(1,1)

    ax.plot(output['time']/365,output['Weight'],'r-')

    ax.set_ylabel('Total Mass Gain (g DW)')

    ax.set_xlabel('Time (years)')

    plt.show()

    return


#%%

if __name__ == "__main__":

    print('Executing my_module.py')

    print('--------------------')



    days, dt, par, InitCond = load_defaults()

    output = run(days, dt, par, InitCond)

    plot_weight(output)

    plot_length(output)

    plot_condition(output)

    plot_masscon(output)



    print('--------------------')




