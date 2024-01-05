TODO:
- Skriva om /protein eller ProteinData.jl, som ny ProteinData.jl
- Funktioner:
    - AA för alla residues (trivialt)
    - ChainID för varje residue (trivialt)
    - Rotations (quat) och locations av residues
    - Residue numbers, baserad på alignment
        - informerad från N-C bond lengths och SEQRES
        - Vi måste få ut SEQRES från pdb-filen också
    

    - PDB utilities för att parsea dataset
    - I/O med H5

    - XYZ i form av Array (3, 3, L) # icke-prioriterat
    - Dihedrals av residues # icke-prioriterat
    - Secondary structure # icke-prioriterat

- Clean dataset:
    - sparar ordnade arrays med sista dimensionen == L, där L är antal residues
    - exempelvis record["locations"] == (3, L)


Model (användning av frames):
Input: geometrin av N residues (representerat på godtyckligt sätt)
Output: 1. Ge positionen av residue N+1 2. Givet Input och positionen, ge vectorpart (rotation) av residue N+1

------

Representation av rotation:
Input: Quaternion som representerar rotationen
1. Skalar quaternionen så att scalarpart == 1 gäller
2. Returnera vectorpart

------

METOD 1:
input: koordinater till en residue: N, CA, C
output: rotation och position
konstant: N-CA, CA-C, <_NCAC

1. position = centroiden av triangeln N, CA, C
2. translaterar triangeln så att den hamnar med centroiden på origo
3 (skeva delen). Räknar ut rotationen som minimerar avståndet mellan triangeln (flyttad till origo) jämfört med en standardtriangel

------

METOD 2:
input: koordinater til alla residues

1. Idealisera koordinaterna, givet bond lengths 
2. Utför METOD 1 (då blir steg 3. mindre skevt)

function metod2(coordinates, bond_lengths):
    coordinates = backbone(idealize_lengths(coordinates, bond_lengths))
    return metod1(coordinates)
end