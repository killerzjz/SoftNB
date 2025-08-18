# SoftNB
SDR-based NB-IoT PHY Signal Processing.

This project offers a prototype of signal processing chain for NB-IoT PHY signal with some sample dataset collected by SDR devices (e.g., USRP, HackRF, RTL Dongle) for test. 

For DL_processing, you can start with running 'DL_frames_processing_chain.py', which contains synchronization, CFO estimation, channel equalization and demodulation. The implementation of these modules is in 'utils.py'

For UL processing, we offer a sample data (TX:HackRF, RX:USRP N210, NRep = 4), you we can start with running 'test_with_NRep4.m', module functions are set in each '.m' file.
