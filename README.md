# Cranfield-LEO-PNT
Repository for Cranfield's ASE 2025-2026 LEO-PNT Group project. Aimed at designing a full post-process receiver chain demonstrator for precise Positioning, Navigation and Timing with (future) real-time capabilities. The initial focus is the development of the core software architecture that enables this functionality.

## Algorithm Structure
The system is structured around a processing algorithm composed of three main stages, each of which is critical to achieving high-accuracy PNT performance.

### 1. Acquistion
This initial part takes care of the signal post processed by Software-Designed Radio (SDR). By means of filters applications, a lower noise signal should be obtained. From there the programmer should be able to extract key (non-encrypted) information, such as time of transmission or relevant satellite information (e.g: ID, for later use to obtain TLEs). Alongside, the code shall be able to identify parts of the signal, highlighting the peaks of synchronisation, which are periodic peaks that are repeated every x milliseconds, which are essential for determining time delays caused by Doppler shifts.

### 2. Tracking
In this stage, the system accurately characterises and tracks the Doppler shift from a noisy signal input. Kalman filtering techniques are applied to smooth the data and enable reliable extrapolation of the satellite’s frequency variations over time.

### 3. Navigation
This final stage implements navigation algorithms derived from established literature. The results obtained here represent the ultimate goal of the project: determining the position of the receiver based on the processed and tracked satellite signal data from LEO.

---
## Rules
### 0. General
- **0.1** All team members shall follow this style guide for all MATLAB code in the project.

### 1. Code Layout

#### 1.1 Functions
- **1.1.1** Use **functions** instead of scripts whenever possible to enhance scalability.  
- **1.1.2** Avoid creating unused or unnecessary functions to keep the codebase simple.

#### 1.2 Indentation
- **1.2.1** Use **4 spaces** per indentation level.  
- **1.2.2** All code blocks (e.g., `if`, `for`, `function`) must be indented consistently.

#### 1.3 Line Length
- **1.3.1** Maximum line length: **75 characters**.  
- **1.3.2** If a line exceeds this limit, split it into multiple lines using the continuation symbol `...` at logical breakpoints.

#### 1.4 Whitespace
- **1.4.1** Avoid **extraneous whitespace**. Follow good practices as shown in MATLAB’s official documentation.

### 2. File Header / Template (Mandatory)

- **2.1** Every MATLAB file (`.m`) shall begin with the standard header from  
  `Algorithm/codingTemplate.m`.  
  - **2.1.1** Copy the header directly from `codingTemplate.m` and ensure all fields are completed or updated appropriately.

### 3. Comments

- **3.1** Comments should be **complete sentences**. Capitalise the first word unless it is an identifier starting with a lowercase letter.  
- **3.2** Begin all comments with `%` followed by a **single space**.  
- **3.3** Separate paragraphs inside a block comment with a line containing a single `#`.  
- **3.4** Each function must include **descriptive comments** explaining its purpose, inputs, and outputs to improve readability.

### 4. Naming Conventions

- **4.1** Use **lowerCamelCase** for variable and function names (e.g., `numIterations`, `sensorData`).  
- **4.2** Boolean variables should clearly indicate true/false states (e.g., `isValid`, `hasConverged`).  
- **4.3** Avoid unclear or misleading names such as:
  - Single-letter names (e.g., `i`, `j`, `k`)  
  - Built-in MATLAB function names (e.g., `disp`, `sum`)

---
## 🔧 Setup
### Clone the repository
```bash
git clone https://github.com/Pablo-TR/Cranfield-LEO-PNT.git
```
Create a folder in the same directory as the *README.md* called **Signals**. This folder will store all sample signal files. It is intentionally excluded from the GitHub repository (see .gitignore) due to licensing and legal restrictions.
