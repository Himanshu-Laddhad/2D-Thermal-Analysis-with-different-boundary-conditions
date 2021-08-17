<h1> 2D-Thermal-Analysis-with-different-boundary-conditions</h1>
<h2>Introduction</h2>
  <p>This module gives thermal analysis of a 2D plate in which 4 different boundary conditions can be applied on each boundary. The boundary conditions available are 
    <ol>
      <li>Dirichilectch Condition</li>
      <li>Homogenous Neuman Condition</li>
      <li>Heterogenous Neuman Condition</li>
      <li>Robbins Condition<\li>
    </ol>
  </p>
<h2>Information Requirement</h2>
<p>
  The basic information that has to be provided to the program include:
  <ul>
    <li>Length</li>
    <li>Width</li>
    <li>Successive ratio in x direction</li>
    <li>Successive ratio in y direction</li>
    <li>Cells in x direction</li>
    <li>Cells in y direction</li>
    <li>Boundary comditions and related data</li>
    <li>Number of iterations</li>
  </ul>
 For graph plotting, the information required is as follows:
  <ul>
    <li>Contour plot is plotted automatically.</li>
    <li>3d Surface plot is plotted automatically.</li>
    <li>Temperature distribution for particular value of x - cell number (in x direction) is required</li>
    <li>Temperature distribution for particular value of y - cell number (in y direction) is required</li>
    <li>Cells in x direction</li>
  </ul>
</p>
<h2>Examples</h2>
  <h3>Meshes</h3>
    <p>Following is a basic uniform square mesh plotted:</p>
      <img src = "https://user-images.githubusercontent.com/63182419/128965926-6bd5b95a-e644-4b33-89ab-32e97afe12f4.png"></img>
<p>Uniform Square Mesh</p>
    <p>Following is a basic non-uniform square mesh plotted:</p>
      <img src = "https://user-images.githubusercontent.com/63182419/128965930-84635c67-54de-4f9d-bf29-3cf066a3917d.png"></img>
<p>Non-Uniform Mesh</p>
  <h3>Graph</h3>
    <p>Following is an example with following conditons:</p>
     <ul>
      <li>Left boundary: Dirichilectch Condition</li>
      <li>Top boundary: Homogenous Neuman Condition</li>
      <li>Bottom boundary: Heterogenous Neuman Condition</li>
      <li>Right boundary: Robbins Condition</li>
     </ul>
   <img src="https://user-images.githubusercontent.com/63182419/129669197-4b47f243-1ded-4f9e-b291-702971015e68.png">
  <break>
  <break>
<footer>
  Countour and 3d surface plot are not available for rectangular meshes
</footer>
