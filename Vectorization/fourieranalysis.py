# Import the required packages
import numpy as np
from scipy.fft import fft, rfft
from scipy.fft import fftfreq, rfftfreq
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import matplotlib.pyplot as plt

# Building a class Fourier for better use of Fourier Analysis.
class Fourier:
  """
  Apply the Discrete Fourier Transform (DFT) on the signal using the Fast Fourier 
  Transform (FFT) from the scipy package.

  Example:
    fourier = Fourier(signal, sampling_rate=2000.0)
  """

  def __init__(self, signal, sampling_rate):
    """
    Initialize the Fourier class.

    Args:
        signal (np.ndarray): The samples of the signal
        sampling_rate (float): The sampling per second of the signal
    
    Additional parameters,which are required to generate Fourier calculations, are
    calculated and defined to be initialized here too:
        time_step (float): 1.0/sampling_rate
        time_axis (np.ndarray): Generate the time axis from the duration and
                              the time_step of the signal. The time axis is
                              for better representation of the signal.
        duration (float): The duration of the signal in seconds.
        frequencies (numpy.ndarray): The frequency axis to generate the spectrum.
        fourier (numpy.ndarray): The DFT using rfft from the scipy package.
    """
    self.signal = signal
    self.sampling_rate = sampling_rate
    self.time_step = 1.0/self.sampling_rate
    self.duration = len(self.signal)/self.sampling_rate
    self.time_axis = np.arange(0, self.duration, self.time_step)
    self.frequencies = rfftfreq(len(self.signal), d = self.time_step)
    self.fourier = rfft(self.signal)
    
  # Generate the actual amplitudes of the spectrum
  def amplitude(self):
    """
    Method of Fourier

    Returns:
        numpy.ndarray of the actual amplitudes of the sinusoids.
    """
    return 2*np.abs(self.fourier)/len(self.signal)

  # Generate the phase information from the output of rfft  
  def phase(self, degree = False):
    """
    Method of Fourier

    Args:
        degree: To choose the type of phase representation (Radian, Degree).
                By default, it's in radian. 

    Returns:
        numpy.ndarray of the phase information of the Fourier output.
    """
    return np.angle(self.fourier, deg = degree)

  # Plot the spectrum
  def plot_spectrum(self, interactive=False):
    """
    Plot the Spectrum (Frequency Domain) of the signal either using the matplotlib
    package, or plot it interactive using the plotly package.

    Args:
        interactive: To choose if you want the plot interactive (True), or not
        (False). The default is the spectrum non-interactive.

    Retruns:
        A plot of the spectrum.
    """
    # When the argument interactive is set to True:
    if interactive:
      self.trace = go.Line(x=self.frequencies, y=self.amplitude())
      self.data = [self.trace]
      self.layout = go.Layout(title=dict(text='Spectrum',
                                         x=0.5,
                                         xanchor='center',
                                         yanchor='top',
                                         font=dict(size=25, family='Arial, bold')),
                              xaxis=dict(title='Frequency[Hz]'),
                              yaxis=dict(title='Amplitude'))
      self.fig = go.Figure(data=self.data, layout=self.layout)
      return self.fig.show()
    # When the argument interactive is set to False:
    else:
      plt.figure(figsize = (10,6))
      plt.plot(self.frequencies, self.amplitude())
      plt.title('Spectrum')
      plt.ylabel('Amplitude')
      plt.xlabel('Frequency[Hz]')
  
  # Plot the Signal and the Spectrum interactively
  def plot_time_frequency(self, t_ylabel="Amplitude", f_ylabel="Amplitude",
                          t_title="Signal (Time Domain)",
                          f_title="Spectrum (Frequency Domain)"):
    """
    Plot the Signal in Time Domain and Frequency Domain using plotly.

    Args:
        t_ylabel (String): Label of the y-axis in Time-Domain
        f_ylabel (String): Label of the y-axis in Frequency-Domain
        t_title (String): Title of the Time-Domain plot
        f_title (String): Title of the Frequency-Domain plot 

    Returns:
        Two figures: the first is the time-domain, and the second is the
                     frequency-domain.
    """
    # The Signal (Time-Domain)
    self.time_trace = go.Line(x=self.time_axis, y=self.signal)
    self.time_domain = [self.time_trace]
    self.layout = go.Layout(title=dict(text=t_title,
                                       x=0.5,
                                       xanchor='center',
                                       yanchor='top',
                                       font=dict(size=25, family='Arial, bold')),
                            xaxis=dict(title='Time[sec]'),
                            yaxis=dict(title=t_ylabel),
                            width=1000,
                            height=400)
    fig = go.Figure(data=self.time_domain, layout=self.layout)
    fig.show()
    # The Spectrum (Frequency-Domain)
    self.freq_trace = go.Line(x=self.frequencies, y=self.amplitude())
    self.frequency_domain = [self.freq_trace]
    self.layout = go.Layout(title=dict(text=f_title,
                                       x=0.5,
                                       xanchor='center',
                                       yanchor='top',
                                       font=dict(size=25, family='Arial, bold')),
                            xaxis=dict(title='Frequency[Hz]'),
                            yaxis=dict(title=f_ylabel),
                            width=1000,
                            height=400)
    fig = go.Figure(data=self.frequency_domain, layout=self.layout)
    fig.show()