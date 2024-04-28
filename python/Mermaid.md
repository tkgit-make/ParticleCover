
```mermaid
graph TD;
    classDef default fill:#333333,stroke:#ffffff,stroke-width:2px, color:#ffffff;
    classDef functionCall fill:#226699,stroke:#ffffff,stroke-width:2px, color:#ffffff;
    classDef loop fill:#884411,stroke:#ffffff,stroke-width:2px, color:#ffffff;
    classDef conditional fill:#446622,stroke:#ffffff,stroke-width:2px, color:#ffffff;



    %% Initialization and Top-Level Loop
    A(Start) --> B(Initialize Parameters: apexZ0, stop, ppl, leftRight);
    B --> C{apexZ0 > -env.trapezoid_edges0?};
    C -- Yes --> D[Adjust apexZ0 to env.trapezoid_edges0];
    C -- No --> End(End);

    %% First Level Loop
    D --> E{number of patches > 0};
    E -- Yes --> F(Adjust z_top_max);
    E -- No --> G[Initialize Variables: z_top_min, complementary_apexZ0, first_row_count, c_corner, nPatchesInColumn, projectionOfCcornerToBeam];
    F --> G;
    G --> H{c_corner > -env.trapezoid_edges.num_layers-1. AND nPatchesInColumn < 100000000 AND projectionOfCcornerToBeam < env.beam_axis_lim};

    %% Inner Loop - Processing patches
    H -- Yes --> I(Update nPatchesInColumn);
    H -- No --> C;
    I --> J(Create Patch alignedToLine);
    J --> K{Iterate over superpoints};

    %% Superpoints Iteration
    K -- Yes --> L(Print superpoint details) --> K;
    K -- No --> M{number of patches > 2};
    M -- Yes --> N(Check Repeat Original Condition);
    M -- No --> O(Calculate projectionOfCcornerToBeam);
    N --> O;

    %% Check Patch Conditions
    O --> P{notChoppedPatch AND more conditions};
    P -- Yes --> Q(Create Complementary Patch);
    P -- No --> R{More iterations needed?};
    R -- Yes --> S(Update z_top_max based on conditions) --> H;
    R -- No --> C;

    %% Complementary Patch Creation
    Q --> T(Calculate Complementary Patch Details) --> U(Adjust for White Space and Overlaps);
    U --> V{Horizontal Shifts Needed?};
    V -- Yes --> W(Apply Shifts and Recalculate Shadows);
    V -- No --> X(Finalize Patch Details);
    W --> X;
    X --> Y{Another Iteration?};
    Y -- Yes --> Z(Adjust Parameters for Next Iteration) --> H;
    Y -- No --> C;

    %% Styling
    class D,E,F,G,H,I,J,K,L,M,N,O,P,Q,R,S,T,U,V,W,X,Y,Z loop;
    class C conditional;
    class J,Q,W functionCall;
```

#```mermaid
graph TD;
    classDef default fill:#f9f,stroke:#333,stroke-width:2px;
    classDef functionCall fill:#ff9,stroke:#333,stroke-width:2px;
    
    %% Initialization 
    A(Start) --> B(Initialize Parameters: apexZ0, stop, ppl, leftRight);
    B --> C[Adjust apexZ0 to env.trapezoid_edges0];
    C --> D{While apexZ0 > -env.trapezoid_edges0?};
    D -->|Yes| E(Initialize Variables: z_top_min, complementary_apexZ0, first_row_count, c_corner);
    E --> F[Calculate z_top_max];
    F --> G{Number of patches > 0?};
    G -->|Yes| H(Adjust z_top_max);
    G -->|No| I;
    H --> I(Inner While Loop: Check conditions);
    I --> J(c_corner > -env.trapezoid_edges.num_layers-1.? AND nPatchesInColumn < Limit AND projectionOfCcornerToBeam < env.beam_axis_lim);
    J -->|Yes| K(nPatchesInColumn++);
    K --> L(Create Patch alignedToLine);
    L --> M{Iterate over superpoints};
    M -->|For each superpoint| N(Print superpoint details);
    N --> O{More superpoints?};
    O -->|Yes| M;
    O -->|No| P(Assess Patch Conditions);
    P --> Q{len-patches > 2?};
    Q -->|Yes| R(Check Repeat Original Condition);
    Q -->|No| S(Set seed_apexZ0);
    R --> S;
    S --> T(Calculate projectionOfCcornerToBeam);
    T --> U(Evaluate Square Patch Alternates);
    U --> V{notChoppedPatch AND Conditions?};
    V -->|No| X(Complementary Patch Conditions);
    X -->|Yes| Y{Create Complementary Patch};

    %% complementary patch creation
    Y --> Z(Calculate Complementary Patch Details);
    Z --> A1{Adjust for White Space and Overlaps};
    A1 --> A2{Check Horizontal Shifts Needed};
    A2 -->|Shifts Needed| A3(Apply Horizontal Shifts and Recalculate Shadows);
    A2 -->|No Shifts| A4(Finalize Patch Details);
    A3 --> A4;
    A4 --> A5{Another Iteration?};
    A5 -->|Yes| D;
    A5 -->|No| A6(End);

    class L,Y,A3 functionCall;

```
