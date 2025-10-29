# Cranfield-LEO-PNT
Repository for Cranfield's ASE 2025-2026 LEO-PNT Group project. Aimed at designing a full post-process receiver chain demonstrator for precise Positioning, Navigation and Timing with (future) real-time capabilities. 

## 🔧 Setup

This project is based around an algorithm, with three components that are key to drive the scope of the project.

### 1. Acquistion
This initial part takes care of the signal post processed by Software-Designed Radio (SDR). By means of filters applications, a lower noise signal should be obtained. From there the programmer should be able to extract key (non-encrypted) information, such as time of transmission or relevant satellite information (e.g: ID, for later use to obtain TLEs). Alongside, the code shall be able to identify parts of the signal, highlighting the peaks of synchronisation, which are periodic peaks that are repeated every x milliseconds, which are essential for determining time delays caused by Doppler shifts.

### 2. Tracking
In this stage, the system accurately characterizes and tracks the Doppler shift from a noisy signal input. Kalman filtering techniques are applied to smooth the data and enable reliable extrapolation of the satellite’s frequency variations over time.

### 3. Navigation

## Rules


## Clone the repository
```bash
git clone https://github.com/Pablo-TR/Cranfield-LEO-PNT.git
```
