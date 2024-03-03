import numpy as np
from behaviour import RegionalBehaviourModel
import pathosim as inf
import matplotlib.pyplot as plt
import os
import pandas as pd

#CREATING THE REGIONAL POPULATION SETUP FOR THE SIMULATION BASED ON ENGLAND DATA

# Actual population data for each region
actual_population_sizes = [9217265, 9002488, 7341196, 6236072, 5934037, 5624696, 5502967, 4835928, 2669941]  # Example data in thousands

# Calculate the total actual population
total_actual_population = sum(actual_population_sizes)

# Your desired total population for the simulation
desired_total_population = 5000

# Calculate the scaling factor
scaling_factor = desired_total_population / total_actual_population

# Calculate the population sizes to use in the simulation
actual_population_sizes = np.array(actual_population_sizes)  # Convert to NumPy array
simulation_population_sizes = np.round(actual_population_sizes * scaling_factor).astype(int)

## Calculate the difference between the desired total population and the sum of rounded population sizes
diff = desired_total_population - sum(simulation_population_sizes)

# Adjust one of the regions to make up the difference (e.g., the first region)
simulation_population_sizes[0] += diff

# Update region parameters to include 9 regions
reg_sizes = simulation_population_sizes
reg_names = ['South East', 'London', 'North West', 'East Anglia', 'West Midlands', 
             'South West', 'Yorkshire & Humber', 'East Midlands', 'North East']
reg_params = []
for reg in reg_names:
    cur_params = dict(
        name=reg,
        n=reg_sizes[reg_names.index(reg)],
        com_contacts=11,
        country_location='england'
    )
    reg_params.append(cur_params)

# Probability matrices from csv
finalized_work_mixing_matrix = np.genfromtxt(f"C:\psim\pathosim\September\Data\White Paper Data - RWMM.csv", delimiter=',')
finalized_com_mixing_matrix = np.genfromtxt(f"C:\psim\pathosim\September\Data\White Paper Data - RCMM.csv", delimiter=',')

def array_to_dict(arr, names):
    """
    Converts an array to the corresponding dictionary structure.

    Parameters:
    arr (numpy.ndarray): Square array
    names (list): List of region names

    Returns:
    dict: Dictionary representation
    """
    n = len(names)
    if arr.shape != (n, n):
        raise ValueError("Array shape must match the length of names.")

    result = {}
    for i in range(n):
        # Overall probability of leaving the home region
        leaving = 1 - arr[i, i]

        # Probabilities of going to other regions
        dests = {}
        for j in range(n):
            if i != j:  # Skip the diagonal
                dests[names[j]] = arr[i, j]
        
        # Normalize the probabilities in dests
        total_dests = sum(dests.values())
        for key in dests.keys():
            dests[key] /= total_dests

        result[names[i]] = {'leaving': leaving, 'dests': dests}

    return result

# Example usage:
names = ['South East', 'London', 'North West', 'East Anglia', 'West Midlands', 
         'South West', 'Yorkshire & Humber', 'East Midlands', 'North East']

arr = finalized_work_mixing_matrix

wmparams = array_to_dict(arr, names)

# Create pop object with new parameters
pop = RegionalBehaviourModel(com_mixing=True, work_mixing=True, 
                             params_com_mixing=finalized_com_mixing_matrix,
                             params_work_mixing=wmparams, 
                             all_reg_params=reg_params)

#BASIC SIMULATION PARAMETERS
region_seed_infections = dict(zip(reg_names, [0, 10, 0, 0, 0, 0, 0, 0, 0]))
multi_reg_params = dict(
    reg_params=pop.reg_pars
)
pars = dict(
    pop_size=desired_total_population,
    pop_type='behaviour_module',
    n_days = 100,
    rand_seed = 0,
    start_day = '2020-01-30',
    enable_multiregion = True,
    multiregion = multi_reg_params
)

#DEFINE PARAMETER SWEEPS. TODO: select appropriate betas.

sample_sizes = [0.000125, 0.00125, 0.0125, 0.125]
sample_frequences = [1, 7, 14, 28]
incubation_periods = [1, 4.5, 14, 28]
betas = [0.004, 0.008, 0.016, 0.032]

#FUNCTIONS FOR DIFFERENT SURVEILLANCE ARCHITECTURES. 

#syndromic surveillance
def syndromic_surveillance(pars, sample_size):

    syndromic_pars = pars.copy()
    syndromic_pars['enable_surveillance'] = True
    syndromic_pars['enable_syndromic_testing'] = True
    syndromic_pars['syndromic_test_percent'] = sample_size #percentage of population to test per day. currently set to replicate threatnet (25% of hospital visits, 30% of hospitals)
    syndromic_pars['hospital_capacity_percent'] = 0.0024 #percentage of hospital capacity to trigger syndromic surveillance - based on UK data.
    syndromic_pars['syndromic_days_to_test'] = 7 #number of days in hospital before testing, arbitrarily set
    syndromic_pars['syndromic_time_to_confirmation'] = 7 #time to confirm a positive test, arbitrarily set
    return syndromic_pars

#random surveillance
def random_surveillance(pars, sample_size, sampling_frequency):

    random_pars = pars.copy()
    random_pars['enable_surveillance'] = True
    random_pars['enable_random_testing'] = True
    random_pars['surveillance_test_percent'] = sample_size #ercentage of population to test per day
    random_pars['random_test_frequency'] = sampling_frequency #frequency of random testing
    
    return random_pars

#regional surveillance
def regional_surveillance(pars, sample_size, sampling_frequency):
    
    regional_pars = pars.copy()
    regional_pars['enable_surveillance'] = True
    regional_pars['enable_regional_testing'] = True
    regional_pars['enable_random_testing'] = True #enable random testing within the region
    regional_pars['surveillance_test_percent'] = sample_size #percentage of population to test per day
    regional_pars['random_test_frequency'] = sampling_frequency #frequency of random testing
    regional_pars['test_regions'] = [1] #region to test
    
    return regional_pars

#age-based surveillance
def age_surveillance(pars, sample_size, sampling_frequency):
    
    age_pars = pars.copy()
    age_pars['enable_surveillance'] = True
    age_pars['enable_age_testing'] = True
    age_pars['surveillance_age_lower'] = 80 #upper age bound
    age_pars['surveillance_test_percent'] = sample_size #ercentage of population to test per day
    age_pars['age_test_frequency'] = sampling_frequency #frequency of random testing
    
    return age_pars

#aiport surveillance
def airport_surveillance(pars, sample_size, sampling_frequency):
    
    airport_pars = pars.copy()
    airport_pars['enable_surveillance'] = True
    airport_pars['enable_airport_testing'] = True
    airport_pars['airport_test_percent'] = sample_size #percentage of population to test per day
    airport_pars['airport_test_frequency'] = sampling_frequency #frequency of random testing
    
    return airport_pars

#FUNCTION FOR RUNNING SIMULATIONS

# Initialize a list to hold the results of each simulation
results_list = []
detection_list = []
run_list = []
cost_list = []
infections_at_detection_list = []

# Initialise mean detection, mean runs, and mean costs
mean_results = []
mean_detection = np.nan
mean_runs = np.nan
mean_costs = np.nan
mean_infections_at_detection = np.nan

#Initialise label
label = ''

def run_simulation(pars, sample_size, sampling_frequency, beta, incubation_period, loops = 100):

    # Run 100 simulations
    for i in range(loops):
        cov = inf.SARS_COV_2(region_seed_infections)

        # Configure the incubation period and beta
        cov.configure_transmission(beta=beta)
        dur = dict(dist='lognormal_int', par1=incubation_period, par2=1.5)
        cov.configure_disease_state_durations(exp2inf=dur)

        # Create run simulation
        sim = inf.Sim(pars=pars, pathogens=[cov], verbose=1, rand_seed=np.random.randint(0, 2**32 - 1, dtype=np.int64))
        sim.initialize()

        # Create airport layers:

        # Create the High-Contact Layer
        p1 = []
        p2 = []
        additional_contacts = 2  # Number of additional contacts for high-contact individuals, aribtrarily set

        # Assuming you want the top X% from each region to be high-contact (i.e. visiting hospitals)
        high_contact_pct = 0.15 # calibrated on https://www.gov.uk/government/statistics/public-experiences-of-and-attitudes-towards-air-travel-2014
        high_contact_individuals = {}

        for reg in reg_params:
            n_people = reg['n']
            high_contact_count = int(n_people * high_contact_pct)
            high_contact_individuals[reg['name']] = list(range(high_contact_count))

            for i in high_contact_individuals[reg['name']]:  # Corrected this line
                for _ in range(additional_contacts):
                    p1.append(i)
                    # Choose contacts randomly or use other criteria
                    contact = (i + 1 + _) % n_people  # Example: choosing them sequentially
                    p2.append(contact)

        #identify airport visitors
        sim.airport_visit_array[p1] = 1

        beta = np.ones(len(p1))  # Transmission risk for the high-contact layer
        layer = inf.Layer(p1=p1, p2=p2, beta=beta)  # Assuming 'inf' is your infection model

        # Add the High-Contact Layer to the Simulation
        sim.people.contacts.add_layer(highcontact=layer)  # Assuming 'Sim' is your simulation object
        sim.reset_layer_pars()

        # Run the simulation
        sim.run()

        # Append results to lists
        results_list.append(sim.results[0]['cum_infections'].values)
        run_list.append(sim.total_runs)
        cost_list.append(sim.total_costs)

        # Only append detection results if it exists
        detection_time = sim.confirmed_detection_time
        if detection_time is not None:
            detection_list.append(detection_time)
            infections_at_detection = sim.results[0]['cum_infections'].values[int(detection_time) - 1]
            infections_at_detection_list.append(infections_at_detection)
            print(f'Simulation {i+1} complete with detection on day {detection_time}')
        else:
            print(f'Simulation {i+1} complete with no detection.')

    # Convert the list to a numpy array for easier calculations
    results_list = np.array(results_list)
    detection_list = np.array(detection_list)
    run_list = np.array(run_list)
    cost_list = np.array(cost_list)
    infections_at_detection_list = np.array(infections_at_detection_list)

    #calculate mean runs and mean costs
    mean_runs = np.nanmean(run_list)
    mean_costs = np.nanmean(cost_list)

    # only calculate mean if there are detections
    if len(detection_list) > 0:
        mean_detection = np.nanmean(detection_list)
        mean_infections_at_detection = np.nanmean(infections_at_detection_list)

    #get confidence intervals and generate label
    intervals = get_confidence_intervals(results_list)
    label = generate_label(sample_size, sampling_frequency, incubation_period, beta)

    #calculate mean results
    mean_results = np.mean(results_list, axis=0)

    # After generating the intervals, plot them:
    days = np.arange(0, pars['n_days'] + 1)
    intervals = get_confidence_intervals(results_list)
    plot_confidence_intervals(intervals, days)
    plot_posterior_realizations(results_list, days, mean_detection, mean_costs)

    #Save to CSV
    df = pd.DataFrame({
        'Day': sim.t,
        'Mean_Infections': mean_results
    })

    filename = os.path.join(output_directory, f"Mean_Results_{label}.csv")
    df.to_csv(filename, index=False)

    return {
        'Simulation_Setup': generate_label(incubation_period, beta),
        'Mean_Detection': mean_detection,
        'Mean_Costs': mean_costs,
        'Mean_Runs': mean_runs,
        'Mean_Infections_at_Detection': mean_infections_at_detection
    }

#FUNCTIONS FOR PLOTTING RESULTS

#define output directory
output_directory = "c:\\Outputs\\Simulations\\White Paper Study\\Results"
if not os.path.exists(output_directory):
    os.makedirs(output_directory)

# Function to get confidence intervals from samples
def get_confidence_intervals(samples, upper=95):
    intervals = np.zeros((3, samples.shape[1]))
    intervals[0, :] = np.percentile(samples, (100-upper)/2, axis=0)
    intervals[1, :] = np.percentile(samples, 50, axis=0)
    intervals[2, :] = np.percentile(samples, upper + (100-upper)/2, axis=0)
    return intervals

# Function for plotting confidence intervals
def plot_confidence_intervals(intervals, days):
    # Extracting the lower bound, median, and upper bound
    lower_bound = intervals[0, :]
    median = intervals[1, :]
    upper_bound = intervals[2, :]
    
    plt.figure(figsize=(10, 6))
    
    # Plot the median
    plt.plot(days, median, label='Median', color='blue')
    
    # Fill the area between the lower and upper confidence intervals
    plt.fill_between(days, lower_bound, upper_bound, color='blue', alpha=0.2, label='95% CI')
    
    # Plot the mean detection day if there was a detection
    if len(detection_list) > 0:
        plt.axvline(mean_detection, color='red', linestyle='--', label=f'Mean Detection Day: {int(round(mean_detection))}')
        cost_text = f"Mean Cost of Program: £{mean_costs:,.2f}"
        sims_text = f"Number of Simulations: {len(results_list)}"
        # Choose a suitable position for the annotations at the top left
        plt.annotate(cost_text, xy=(0.025, 0.80), xycoords='axes fraction', fontsize=10, color='green')
        plt.annotate(sims_text, xy=(0.025, 0.75), xycoords='axes fraction', fontsize=10, color='green')

    plt.xlabel('Day')
    plt.ylabel('Number of Infections')
    plt.title('Number of Infections with 95% Confidence Intervals')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plot_filename = os.path.join(output_directory, f"{label}.png")
    plt.savefig(plot_filename, dpi=300)
    plt.close()

    return plot_filename

#Function for plotting posterior realisations
def plot_posterior_realizations(results_list, days, mean_detection=None, mean_costs=None):
    # Convert results_list to a NumPy array for easier manipulation
    results_array = np.array(results_list)
    
    plt.figure(figsize=(10, 6))
    
    # Plot each realization (posterior sample)
    for realization in results_array:
        plt.plot(days, realization, alpha=0.2, color='blue')
    
    # Plot the mean detection day if there was a detection
    if mean_detection is not None:
        plt.axvline(mean_detection, color='red', linestyle='--', label=f'Mean Detection Day: {int(round(mean_detection))}')
        cost_text = f"Mean Cost of Program: £{mean_costs:,.2f}"
        plt.annotate(cost_text, xy=(0.025, 0.80), xycoords='axes fraction', fontsize=10, color='green')
    
    plt.xlabel('Day')
    plt.ylabel('Number of Infections')
    plt.title('Posterior Realizations of Number of Infections')
    plt.grid(True)
    plt.tight_layout()
    
    plot_filename = os.path.join(output_directory, f"Posterior_{label}.png")
    plt.savefig(plot_filename, dpi=300)
    plt.close()

    return plot_filename


# Generate a label for the simulation setup based on the parameters
def generate_label(sample_size, sampling_frequency, incubation_period, beta):
    label = f'Sample size: {sample_size}, Sampling frequency: {sampling_frequency}, Incubation period: {incubation_period}, Beta: {beta}'
    return label

#RUN AND SAVE SIMULATIONS

# Initialize a list to store summary statistics or other results
columns = ['Surveillance_Type', 'Sample_Size', 'Sampling_Frequency', 'Beta', 'Incubation_Period', 'Mean_Detection', 'Mean_Costs', 'Mean_Runs', 'Mean_Infections_at_Detection']
summary_results = pd.DataFrame(columns=columns)

# Define a list of surveillance types
surveillance_types = ['syndromic', 'random', 'age', 'regional', 'airport']

# Loop over all combinations of sample sizes, frequencies, regions, etc.
for surveillance_type in surveillance_types:
    for sample_size in sample_sizes:
        for sampling_frequency in sample_frequences:
            for beta in betas:
                for incubation_period in incubation_periods:

                    # Choose the surveillance architecture based on the current type
                    if surveillance_type == 'syndromic':
                        current_pars = syndromic_surveillance(pars, sample_size, sampling_frequency)
                    elif surveillance_type == 'random':
                        current_pars = random_surveillance(pars, sample_size, sampling_frequency)
                    elif surveillance_type == 'age':
                        current_pars = age_surveillance(pars, sample_size, sampling_frequency)
                    elif surveillance_type == 'regional':
                        current_pars = regional_surveillance(pars, sample_size, sampling_frequency)
                    elif surveillance_type == 'airport':
                        current_pars = airport_surveillance(pars, sample_size, sampling_frequency)

                    # Run the simulation
                    result = run_simulation(current_pars, sample_size, sampling_frequency, beta, incubation_period)

                    new_row = pd.DataFrame({
                        'Surveillance_Type': [surveillance_type],
                        'Sample_Size': [sample_size],
                        'Sampling_Frequency': [sampling_frequency],
                        'Beta': [beta],
                        'Incubation_Period': [incubation_period],
                        'Mean_Detection': [result['Mean_Detection']],
                        'Mean_Costs': [result['Mean_Costs']],
                        'Mean_Runs': [result['Mean_Runs']],
                        'Mean_Infections_at_Detection': [result['Mean_Infections_at_Detection']]
                    })

                    summary_results = pd.concat([summary_results, new_row], ignore_index=True)

# Save the summary results to a CSV file
summary_results.to_csv(os.path.join(output_directory, 'Final_Results.csv'), index=False)
