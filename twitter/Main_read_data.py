import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timedelta, date
from Optimize_SIR2 import SIR_rhs, SIR_Model, error_sum_of_squares, Optimize_SIR
from GroupCount2 import groupcount
import csv
from Datenum import days_from_zero
#import tweetnlp
import os

filename = 'Epidemic_100-101_Python_Fit_Values.csv'
headers = ['Tweet', 'Topic','datetime', 'Rsquare', 'Alpha', 'Beta','1/Alpha','Alpha/Beta','Mu',"Sum of Followers","EpidemicID"]

# Check if file exists and is empty
if not os.path.isfile(filename) or os.stat(filename).st_size == 0:
    write_headers = True
else:
    write_headers = False


for i in range(100, 101):
    trial = str(i)
    epidemic_id = trial
    base_dir = f"epidemics_{trial}/csv/"
    # Read data from CSV files
    raw_data = pd.read_csv(base_dir + "level1.csv")
    raw_data = pd.read_csv(base_dir + 'level1.csv', sep='\t')
    tweet = raw_data['text'].iloc[-1]
    tweet_data = pd.read_csv(base_dir + "original.csv")
    sentiment_data = pd.read_csv(base_dir + "original_sentiment.csv")

    # Extract sentiment
    sentiment = sentiment_data.iloc[:, -3:].values[0]
    # Split sentiment into a list of values - I only need the sentiment values
    # Get the first (and only) element of the array and split it into a list of values
    sentiment_values = sentiment[0].split('\t')

    # Convert the last three values to floats
    S = [float(value) for value in sentiment_values[-3:]]




    # Setup Susceptible (Followers), Infected (Retweets), time classes
    raw_data.columns = raw_data.columns.str.strip()
    time = pd.to_datetime(raw_data['timestamp'], format="%Y-%m-%dT%H:%M:%S.000Z")
    # Writing the Time Array to a CSV file for troubleshooting
    time = time.dt.tz_localize('UTC').dt.tz_convert('America/Denver')
    #time.to_csv('time.csv', index=False)
    followers = raw_data['followers']
    tweet = tweet_data.iloc[:, 5:]


    # Assuming time, tweet, followers are all lists or arrays.
    sir = {} # SIR dictionary stores key parameters 
    sir['total_tweets'] = len(time) #Same as MATLAB
    ''' 
    def seconds_from_zero(dates):
        # Define the reference datetime as January 1, 0001, 00:00:00
        ref_datetime = datetime(1, 1, 1)
        # Convert to datetime if not already        
        # Calculate the difference in seconds
        delta = [d - ref_datetime.total_seconds() for d in dates]
        
        # Return the number of seconds, plus one year to account for the year 0000
        # timedelta.total_seconds() returns the total number of seconds contained in the duration
        return [d.total_seconds() + 31536000 for d in delta]  # 31536000 is the number of seconds in a non-leap year

    def seconds_to_datetime(seconds, timezone_str):
        ref_datetime = datetime(1, 1, 1)  # Reference date
        dt = ref_datetime + timedelta(seconds=seconds-31536000)  # Convert seconds to timedelta, subtract the year 0000
        dt = dt.replace(tzinfo=timezone.utc)  # Set initial timezone to UTC

        # Convert to the desired timezone
        desired_tz = pytz.timezone(timezone_str)
        dt = dt.astimezone(desired_tz)

        return dt
    '''


    sir['start'] = min(time)
    sir['stop'] = max(time)
    print(f"First datetime: {sir['start']}")
    print(f"Last datetime: {sir['stop']}")



    sir['run_time'] = sir['stop'] - sir['start']  # units = seconds. Difference betweem stop and start
    sir['run_time'] = sir['run_time'].total_seconds() / (24 * 60 * 60)
   # --------------------Generating Histogram----------------------------------------
    # Computing histogram with 'Freedman-Diaconis' rule
    b,edges = np.histogram(followers, bins='fd')
    # Prints out number of bins, should be 111. 
    print(len(b))
    # ----------------------------------------------------------------------
    #sir['timestep'] = len(b)
    sir['timestep'] = 48 #Stopgap measure. Fd likely not good use in Python, other methods (Scott) align much better with MATLAB
    sir['del_t'] = sir['run_time'] / sir['timestep']

    binned = {}
    binned['time'], binned['Susc'], binned['Inf'] = groupcount(
        time, tweet, followers, sir['timestep'], sir['del_t'], sir['start'], sir['stop'], sir['total_tweets'])

    #Deleting unnecessary variables

    # ------------------------------------Performing SIR Optimization--------------------------------------
    # Run the SIR model with the Data
    SIR, ESS, SIR_values, SIR_time, Rsquared = Optimize_SIR(binned['time'], binned['Inf'], binned['Susc'], sir['run_time'])
    T = [SIR_values[3], SIR_values[0], SIR_values[1], 1/SIR_values[1], SIR_values[0]/SIR_values[1], SIR_values[2], np.sum(followers)] 
    T = [x.item() if isinstance(x, np.ndarray) else x for x in T]
    #----------------------------------------Topic, Tweet, Datetime------------------------------------------
    '''''
    tweet = raw_data['text'].iloc[-1]
    model = tweetnlp.load_model('topic_classification', multi_label=True)
    nlp_topic = model.topic(tweet)
    datetime_tweet = min(time) #grabbing last datetime 
    data = []
    data.append(tweet)
    data.append(nlp_topic)
    data.append(datetime_tweet)
    data.extend(T) 
    data.append(epidemic_id)
    '''
    # --------------------------------------Plotting SIR Model & Histogram---------------------------------------------------------
    # Plotting - worry about later 
    ''''
    plt.figure()
    plt.hist(binned['time'], binned['Inf'], bins = len((b)), facecolor='y')
    plt.plot(SIR_time, SIR[:,1], 'r', linewidth=1.5)
    plt.legend(['Infected Data', 'Infected SIR Model'], loc='upper right')
    plt.xlabel('Days')
    plt.ylabel('Number of People')
    plt.show()
    '''

    # Assuming binned_time and binned_inf are the equivalent data in Python
    plt.bar(binned['time'],binned['Inf'], facecolor='y')
    # Histogram without SIR model  
    plt.legend(['Infected Data', 'Infected SIR Model'], loc='upper right')
    plt.xlabel('Days')
    plt.ylabel('Number of People')
    plt.show()
    # SIR Model Graph
    plt.plot(SIR_time, SIR[:, 1], 'r', linewidth=1.5)
    plt.xlabel('Days')
    plt.ylabel('Number of People')
    plt.legend(['Infected SIR Model'])
    plt.show() 
    '''
    # --------------------------------------------CSV FILE GENERATION---------------------------------------------
    with open(filename, 'a', newline='') as f:
        writer = csv.writer(f)
        if write_headers:
            writer.writerow(headers)
            write_headers = False
        writer.writerow(data)
    #plt.show()
    '''

#---------------------------------MATLAB Comparison Code-------------------------------------------
'''
%% This script will read the data into a file and sort it to prepare for SIR Model

% Nick Duncan August 2021
clear
clc
%% Read data file 
%% Reading in two separate files:level1.csv and original.csv
for i = 100:100
trial       = i;
begin       = 'epidemics_';
trial       = num2str(trial);
ending      = '/csv/level1.csv';
file        = append(begin,trial,ending);
data.opts   = detectImportOptions(file);
data.raw    = readtable(file,data.opts);
ending      = '/csv/original.csv';
file        = append(begin,trial,ending);
text.opts   = detectImportOptions(file);
text.tweet  =  readtable(file,text.opts);
%% Read in the Sentiment file of the original message
ending              = '/csv/original_sentiment.csv';
file                = append(begin,trial,ending);
data.opts           = detectImportOptions(file);
data.score          = readtable(file,data.opts);
numColumns = width(data.score); % get the number of columns
sentiment = table2array(data.score(1, numColumns-2:numColumns));
% sentiment           = table2array([data.score(1,8),data.score(1,9),data.score(1,10)]);
clear begin trial ending file
%% Setup Susceptible (Followers), Infected (Retweets), time classes
time        = table2cell(data.raw(:,2));
time        = datetime(time,'TimeZone','America/Denver','Inputformat','yyyy-MM-dd''T''HH:mm:ss.SSSZ');
followers   = table2array(data.raw(:,4));
tweet       = text.tweet(:,6:end);
clear text data
%% divide up the results for the SIR model
sir.total_tweets = length(time);
sir.start        = datenum(time(end));
sir.stop         = datenum(time(1));
sir.run_time     = sir.stop - sir.start; % units = days
[b, edges]   = histcounts(followers,'BinMethod','fd');
sir.timestep = length(b);
sir.del_t = (sir.run_time/sir.timestep); 
[binned.time binned.Susc binned.Inf] = groupcount(time, tweet, followers,...
    sir.timestep,sir.del_t,sir.start,sir.stop,sir.total_tweets);
clear b edges
% figure();
% bar(binned.time,binned.Inf,'facecolor','y');
% hold on
% plot(binned.time,binned.Susc,'-r','linewidth',1.5);
%% Run the SIR model with the Data
[SIR,ESS,SIR_values,SIR_time,Rsquared] = Optimize_SIR(binned.time, binned.Inf, binned.Susc,sir.run_time);
figure
bar(binned.time,binned.Inf,'facecolor','y');
hold on
plot(SIR_time,SIR(:,2),'r','linewidth',1.5);
legend('Infected Data', 'Infected SIR Model','location','Northeast');
xlabel('Days');
ylabel('Number of People');
%% Copy metrics to a file
T = [SIR_values(4) SIR_values(1) SIR_values(2) 1/SIR_values(2) SIR_values(1)/SIR_values(2) SIR_values(3) sum(followers)]; 
S = [sentiment(1) sentiment(2) sentiment(3)];
% Concatenate T and S into a single row vector
combined = [T, S];

filename = 'SIR_Sentiment_100_100Runs.csv';

%% Write the combined SIR & Sentiment vector to a csv file
writematrix(combined, filename, 'WriteMode', 'append');
%% Separate SIR Value File
% filename1 = 'politics_SIR_values.csv';
% writematrix(T,filename1,'WriteMode','append');
%% Tweet File
filename2 = 'politics_tweets_100_100Runs.csv';
writetable(tweet,filename2,'WriteMode','append', 'WriteVariableNames',false);
%% Separate Sentiment File
% filename3 = 'politics_sentiment.csv';
% writematrix(S,filename3,'WriteMode','append');
end
fclose('all');
'''