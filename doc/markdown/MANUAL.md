\section Contents
<ul>
	<li>
		\ref sect_fluid_eq
		<ul>
			<li>
				 \ref sect_physical
				<ul>
					<li> \ref sect_ddr </li>
					<li> \ref sect_elecpotfield </li>
					<li> \ref sect_rescaling </li>
				</ul>
			</li>
			<li> \ref sect_bic </li>
			<li>
				 \ref sect_numeric
				<ul>
					<li> \ref sect_dde </li>
					<li> \ref sect_dpe </li>
				</ul>
			</li>
			<li> \ref sect_timestepping </li>
		</ul>
	</li>
	<li>
		\ref sect_refinement
		<ul>
			<li> \ref sect_overview </li>
			<li> \ref sect_size_refinement </li>
			<li> \ref sect_criterion </li>
			<li> \ref sect_curvature </li>
			<li> \ref sect_fishpack </li>
			<li> \ref sect_Conclusions </li>
		</ul>
	</li>
	<li>
		\ref sect_software
		<ul>
			<li> \ref sect_ARCoS_overview </li>
			<li>
				 \ref sect_IO
				<ul>
					<li> \ref sect_sim </li>
					<li> \ref sect_input </li>
					<li> \ref sect_parameter </li>
					<li> \ref sect_needle </li>
					<li> \ref sect_output </li>
					<li> \ref sect_structure </li>
					<li> \ref sect_source </li>
				</ul>
			</li>
		</ul>
	</li>
</ul>

\section sect_fluid_eq Fluid model

\subsection sect_physical Physical model

\subsubsection sect_ddr  Drift-diffusion-reaction equations

[1]: http://homepages.cwi.nl/~ebert/CaroJCP06.pdf		"Montijn et al."
[2]: http://homepages.cwi.nl/~ebert/LuqueAPL07.pdf		"Luque et al."
[3]: http://homepages.cwi.nl/~ebert/StrBranchLuque2011.pdf	"Luque et al.(2011)"
[4]: http://homepages.cwi.nl/~ebert/JCP-Li-12.pdf 		"Li et al."
[5]: http://arxiv.org/abs/1301.1552				"Teunissen et al.".
[6]: http://							"PhD Thesis Wormeester(to appear in 2013)"
[7]: http://alexandria.tue.nl/repository/books/598717.pdf	"PhD Thesis Montijn"
[8]: http://homepages.cwi.nl/~willem/				"Monograph by Hundsdorfer and Verwer"
[9]: http://www2.cisl.ucar.edu/resources/legacy/fishpack	"Fishpack"
[10]: http://www2.cisl.ucar.edu/resources/legacy/fishpack90	"Fish90"
[11]: http://www.amazon.com/books/dp/3540194622			"Y. P. Raizer, Gas Discharge Physics"
[12]: http://www.stack.nl/~dimitri/doxygen/manual/index.html	"Doxygen"


In a fluid model of a streamer, we replace the individual particles in the
system by a density function \f${n}(\mathbf{r},t)\f$. 
The temporal evolution of this density function is governed by the physical
processes of the system and this model takes the form of a set of partial
differential equations (PDEs).
The derivation of this so-called classical streamer model starts from the
<a name="continuity"><i>continuity</i></a> equation.  For particle species \f$i\f$, we have:
\f[
 \frac{\partial n_i(\mathbf{r},t)}{\partial t} + \nabla \cdot \mathbf{j}_i(\mathbf{r},t) = S_i(\mathbf{r},t).
\f]

Here \f$S_i(\mathbf{r},t)\f$ represents the total of all sources and sinks of
species \f$i\f$.  \f$\mathbf{j}_i(\mathbf{r},t)\f$ is the term for the particle
current density of species \f$i\f$.
Particles can drift and diffuse as described by the following
<a name="current_dens"><i>expression</i></a> for the particle current density
\f$\mathbf{j}_i\f$:
\f[
 \mathbf{j}_i(\mathbf{r},t) = \mu_i n_i(\mathbf{r},t) \mathbf{E}(\mathbf{r},t) - D_i \nabla n_i (\mathbf{r},t).

\f]
In the above equation, the first term represents the particle drift due
to the electric field, with \f$\mu_i\f$ the mobility coefficient of species \f$i\f$.
The second term represents the diffusion of particles due to the spatial gradient in particle densities with diffusion coefficient \f$D_i\f$.  These equations can be derived from the Boltzmann equation

On the timescales involved, we consider only electrons to be mobile, while ions and neutrals remain stationary, which means that for heavy species, the [continuity equation](#continuity) is reduced to
\f[
  \frac{\partial n_i(\mathbf{r},t)}{\partial t} = S_i(\mathbf{r},t).
\f]

The sources and sinks in the first equation play a very important role in the dynamics of the streamer.  In this model, the sources and sinks correspond to reactions between the different charged and neutral species present in the gas.  The source term due to a single reaction is the product of the densities of the species involved in the reaction and the field-dependent rate coefficient for that reaction.

As an example, the <a name="reaction_imp_ion"><i>impact ionization reaction</i></a>
\f[
 e^- + \mathbf{N}_2 \rightarrow 2 e^- + \mathbf{N}_2^+
\f]
is modeled by
\f[
 S_{ionization} = k_{ion}(|\mathbf{E}|) n_e [\mathbf{N}_2].
\f]
Here \f$n_e\f$ is the local electron density, \f$[\mbox{N}_2]\f$ the density of
\f$\mbox{N}_2\f$ and \f$k_{ion}\f$ the reaction coefficient for impact ionization
depending on the magnitude of the local electric field.
The value of \f$k\f$ can be determined in different ways, from experiments,
theoretical calculations or simulations.
The traditional approximation suggested by <a name="townsend_ionization"><i>Townsend</i> </a>
uses an empirical expression for the impact ionization term, see the 
[book by  Y. P. Raizer, "Gas Discharge Physics"][11]:
\f[
 \frac{dn_e}{dt} = n_e \mu_e |\mathbf{E}| \alpha_0 e^{-E_0 / |\mathbf{E}|},
\f]
where \f$\mu_e\f$ is the electron mobility coefficient,
\f$\mathbf{E}\f$ is the local electric field and \f$\alpha_0\f$ and \f$E_0\f$ are
parameters that can be determined by fitting experimental data.
In gases that contain an electronegative admixture, such as \f$\mathbf{O}_2\f$, the process of
attachment can provide a sink for the electron density through the following
reactions:
\f[
 e^- + \mathbf{O}_2 \rightarrow \mathbf{O} + \mathbf{O}^-
\f]
\f[
 e^- + \mathbf{O}_2 + \mathbf{O}_2 \rightarrow \mathbf{O}_2^- + \mathbf{O}_2
\f]
The first attachment process is dissociative attachment, the second an example
of a 3-body attachment (a 3-body attachment can also occur with an oxygen and nitrogen
molecule).
In the case of the 3-body attachment, the reaction rate scales with the square of the oxygen density:
\f[
 S_{3-body-att} = k_{3-body-att}(|\mathbf{E}|) n_e [\mathbf{O}_2]^2.
\f]
Further ionization losses can occur via one or more recombination processes,
but these typically have a timescale that is much longer than the timescale of
streamer development and propagation and are therefore primarily interesting for
the evolution of the charge density after a streamer discharge,
as discussed in [PhD Thesis Wormeester(to appear in 2013)][6], Chapter 5.



In gases with attachment, detachment may occur, resulting in an additional
source of electrons.
In gases that contain both nitrogen and oxygen, the photoionization process
provides a non-local source of electrons.
Since photoionization is non-local, it can not be modelled by simple reaction
equations such as the ones for impact ionization.
Instead, the local contribution of photoionization is calculated by spatially
integrating contributions from the entire domain.
The commonly used model for photoionization and the approximations made to make
this model suitable for simulation are discussed in
[PhD Thesis Wormeester(to appear in 2013)][6], Chapter 3, section Photoionization.

The reaction model for streamer simulations can be very minimal or very extended,
with many species and reactions, including metastables and various excited states.
The complexity of the reaction model depends on the purpose of the simulations.
For negative streamers in nitrogen, a model containing no more than three species
(\f$\mbox{e}^-\f$, \f$\mbox{N}_2\f$ and \f$\mbox{N}_2^+\f$) and one reaction
(impact ionization)
is sufficient to simulate the dynamics of the streamer head, see [Montijn et al][1].
For more detailed studies of the streamer chemistry, the reaction model should be
as complete as possible.


\subsubsection sect_elecpotfield Electric potential and field

The streamer evolves under the influence of an electric field, which consists
of an externally applied electric field and the electric field generated by
space charges.
These space charges are present at the head of the streamer as well as on the edge of the streamer channel.
For the further propagation of the streamer, the enhanced electric field in front of the streamer, generated by the space charge in the streamer head is essential.
We compute the net <a name="charge_dens"><i>charge density</i></a>
\f$q(\mathbf{r},t)\f$:
\f[
 q(\mathbf{r},t) = e \sum_i q_i n_i(\mathbf{r},t),
\f]
where for species \f$i\f$, \f$n_i\f$
denotes the density function of these species and \f$q_i\f$ the charge of
a particle in units of the electron charge \f$e\f$.
From this we compute the potential by solving the <a name="Poisson"><i>Poisson</i> equation
\f[
 \nabla^2 \phi(\mathbf{r},t) = \frac{q(\mathbf{r},t)}{\epsilon_0}
\f]
and the electric field
\f[
 \mathbf{E}(\mathbf{r},t) = -\nabla \phi(\mathbf{r},t).
\f]

\subsubsection sect_rescaling Rescaling to dimensionless units

The classical fluid model for streamers can be rescaled to dimensionless units and
it is with these units that the code used in this documentation works.
We refer the interested reader to the
[PhD thesis of Gideon Wormeester (to appear in 2013)][6].
From the [Townsend](#townsend_ionization) approximation for ionization,
a characteristic field and length scale
emerges: \f$E_0\f$ and \f$l_0 = \alpha_0^{-1}\f$, respectively.
The characteristic velocity follows from the drift velocity of electrons
in the characteristic field,
\f[
  E_0: v_0 = \mu_e E_0.
\f]

The characteristic number density follows from the [Poisson](#Poisson) equation. 
Values for \f$\alpha_0\f$, \f$E_0\f$ and \f$\mu_e\f$ were obtained from
[PhD Thesis Montijn][7]
and are at standard temperature and pressure:
\f{eqnarray*}
 \alpha_0 & \simeq &  4332          \quad      \mbox{cm}^{-1}\\
 E_0      & \simeq &  2 \times 10^5 \quad \mbox{V} \mbox{cm}^{-1}\\
 \mu_e    & \simeq &  380           \quad \mbox{~cm}^2 \mbox{V}^{-1} \mbox{s}^{-1}.
\f}

When we insert these values in the characteristic scales, we obtain the values with which to rescale the equations:
\f{eqnarray*}
 l_0 & \simeq & 2.3 \times 10^{-4} \quad \mbox{~cm}\\
 t_0 & \simeq & 3.0 \times 10^{-12} \quad \mbox{~s}\\
 n_0 & \simeq & 4.7 \times 10^{14} \quad \mbox{~cm}^{-3}\\
 D_0 & \simeq & 1.8 \times 10^{4} \quad \mbox{~cm}^2 \mbox{s}^{-1}.
\f}
We can now make the appropriate substitutions
(\f$t^d = t / t_0\f$ and similarly for the other variables;
the superscript \f${\;}^d\f$ will be used to indicate that a variable is
in dimensionless form, where this is not clear from the context.
For clarity of reading, the superscript \f${\;}^d\f$ will be omitted where it is clear that variables are
dimensionless) to obtain the classical fluid equations in
<a name="continuity_dimless"><i>dimensionless continuity</i></a> form:
\f[
 \partial_{t} + \nabla \cdot \mathbf{j}_i = S_i,
\f]
where \f$t\f$ is the dimensionless time, \f$\mathbf{j}_i\f$ the dimensionless
particle density current for species \f$i\f$ and \f$S_i\f$ the dimensionless source term
for species \f$i\f$.
\f$S_i\f$ is obtained by rewriting reaction equations such as the
[impact ionization reaction](#reaction_imp_ion) equation in dimensionless form,
where we remark that all rate-coefficients should also be rescaled.
The particle density current \f$\mathbf{j}_i\f$ is obtained by rescaling
[expression](#current_dens) into
<a name="current_dens_dimless"><i>equation</i></a>
\f[
 \mathbf{j}_i = -\mu_i n_i \mathbf{E} - D_i \nabla n_i,
\f]
where \f$\mathbf{E}\f$ is the dimensionless electric field and \f$n_i\f$,
\f$D_i\f$ and \f$\mu_i\f$ are the dimensionless particle density,
diffusion coefficient and mobility respectively of species \f$i\f$.
We find that in dimensionless units \f$\mu_i\f$ is equal to 1 while for heavy particles
\f$\mu_i\f$ is taken as 0, since heavy particles are assumed to be stationary in
this model.
The [dimensionless current density equation](#current_dens_dimless) can therefore be simplified to
\f[
 \mathbf{j}_e = -n_e \mathbf{E} - D_e \nabla n_e
\f]
for electrons and \f$ \mathbf{j}_i = 0 \f$ for heavy particles.
The expression for the [charge density equation](#charge_dens) \f$q\f$, is rescaled to
\f[
 q(\mathbf{r},t) = \sum_i q_i n_i(\mathbf{r},t).
\f]
The [Poisson](#Poisson) equation is rescaled to
\f[
 \nabla^2 \phi = q.
\f]
We remark that although the code
described here internally
works with the dimensionless equations and variables described in this section,
all results are presented in regular units unless otherwise noted.
Input parameters for the simulation code are expected to be in dimensionless units.
Finally we note that the rescaling to dimensionless units does not change
the structure of the equations, it is merely a rescaling to a different
set of units, where the dimensionless units yield a set of equations where some
constants (such as \f$e\f$, \f$\epsilon_0\f$, \f$\mu_e\f$) become unity.

\subsection sect_bic Boundary and initial conditions


We consider a cylindrical computational domain with coordinates:
\f[
(r,z,\theta) \in (0,L_r) \times (0,L_z) \times (0,2\pi).
\f]
Although the code described here is capable of performing full 3D calculations,
we assume cylindrical symmetry to greatly simplify the computations.
For any spatially dependent function \f$f(r,z,\theta)\f$, we assume:
\f$\partial_{\theta} f(r,z,\theta) = 0\f$.
Consequently, the coordinate system for our computations is limited to
\f$(0,L_r) \times (0,L_z)\f$.
We consider a setup with a powered electrode at \f$z = L_z\f$ and a grounded
electrode at \f$z = 0\f$.
If the powered electrode is a plate, the following boundary conditions are used for
the electric potential \f$\phi(r,z,t)\f$:

\f[
 \begin{array}{llll}
 \forall z \; & \partial_r \phi(0,z,t) & = & 0\\
 \forall r \; & \phi(L_r,z,t) & = & 0\\
 \forall z \; & \phi(r,0,t)   & = & 0\\
 \forall r \; & \phi(r,L_z,t) & = & \phi_0 
 \end{array}
\f]

with \f$\phi_0\f$ the potential applied to the powered electrode.
If the powered electrode is a needle protruding from a plate, the needle has the
same potential \f$\phi_0\f$ as the plate.

<img src="http://md-wiki.project.cwi.nl/images/figure2.png" width="240px" />
\par Figure 1. Schematic of the computational setup.

In <a name="setup">Figure 1.</a>, the shaded rectangle represents the computational domain for
the fluid equations, the thick horizontal lines the two planar electrodes
with the needle and its parameters depicted at the anode.
The area between the two planar electrodes is the computational domain
for the Poisson equation.
The needle is simulated by a single point charge, \f$Q\f$, chosen such
that \f$\phi =\phi_0\f$ in the point \f$P\f$, which is the tip of the needle.
The calculation assumes cylindrical symmetry around the needle axis
represented by the dashed-dotted line.

For the density equations, we use homogeneous Neumann conditions on all edges:
\f[
 \partial_r n(0,z,t) = \partial_r n(L_r,z,t) = \partial_z n(r,0,t) = \partial_z n(r,L_z,t) = 0,
\f]
where we remark that if the powered electrode is a needle, the computational domain
for the density equations is smaller than the computational domain for the Poisson
equation and the \f$L_z\f$ values for both domains are not equal.
This difference is a requirement of the numerical implementation of the needle
electrode and is further detailed in section \ref sect_needle.

While the boundary conditions mentioned above are the ones used by Wormeester,
the code that was used can also handle different choices of boundary conditions:
both homogeneous Neumann and homogeneous Dirichlet boundary conditions are available
for the top (\f$z = L_z\f$), bottom (\f$z = 0\f$) and right (\f$r = L_r\f$)
edges of the domain for both the densities and the potential.
The Neumann condition on the central axis of the cylindrical domain is required
for symmetry reasons.

As initial conditions for particle densities, two types of seeds are implemented in
the code. A homogeneous seed, with a constant density over the entire domain and a
Gaussian seed of the form
\f[
 n(r,z,0) = n_{max} \mbox{exp}(-\frac{r^2 + (z - z_0)^2}{\sigma^2}).
\f]
Here \f$z_0\f$ specifies the \f$z\f$-coordinate of the maximum of the seed
(which is located on the symmetry axis with \f$r = 0\f$), where the density is
\f$n_{max}\f$. \f$\sigma\f$ is a measure of the radius of the seed, it is the distance
at which the density drops to \f$e^{-1}\f$ of the maximum value.

In typical streamer simulations, a seed of electrons and positive ions is placed at the tip of the needle to initiate the discharge. Other than these Gaussian seeds and the neutral background gas, initial particle densities are zero with the possible exception of added background ionization, a homogeneous density of negative and positive ions. The initial distribution of electrons and ions is charge neutral at every point of the domain.

\subsection sect_numeric Numerical method

The physical equations in section \ref sect_fluid_eq are to be solved numerically.
The computational code we have applied for this uses finite volume methods to solve
a discretized version of the physical equations.
Here we give a basic summary of the numerical technique used.
For more details, the reader is referred to the work of [Montijn et al.][1],
upon which the current code is based.

\subsubsection sect_dde Discretization of density equations

The [dimensionless continuity ](#continuity_dimless) and the [dimensionless current density](#current_dens_dimless) equations
are discretized using finite volume methods and solved on a uniform rectangular grid with cells:
\f[
 C_{ij} = [(i - 1) \Delta r, i \Delta r] \times [(j - 1) \Delta z, j \Delta z]\left(i = 1 , \cdots , \frac{L_r}{\Delta r}, j = 1 , \cdots , \frac{L_z}{\Delta z}\right),
\f]
where \f$L_r\f$ and \f$L_z\f$ are the \f$r\f$- and \f$z\f$-dimensions of the grid and \f$\Delta r\f$ and \f$\Delta z\f$ the size of a cell in \f$r\f$- and \f$z\f$-direction, respectively. Particle density distributions are represented by their value in the cell center, which can be seen as an average over the cell. For some species \f$n\f$, we use \f$n_{i,j}\f$ to denote the density at the cell center \f$C_{ij}\f$. For sake of clarity of notation we omit the superscript \f$^d\f$ indicating that variables are in dimensionless units.

The discretized continuity equations in cylindrical coordinates, with cylindrical symmetry (\f$\partial_{\theta} f = 0\f$) assumed, have the following form:
\f[
 \begin{array}{ll}
 \frac{d n_{i,j}}{d t} = & \frac{1}{r_i \Delta r} \Big(r_{i - \frac{1}{2}} F^a_{i - \frac{1}{2},j} - r_{i + \frac{1}{2}} F^a_{i + \frac{1}{2},j} + r_{i - \frac{1}{2}} F^d_{i - \frac{1}{2},j} - r_{i + \frac{1}{2}} F^d_{i + \frac{1}{2},j}\Big) + \\
                         & \frac{1}{\Delta z} \Big(F^a_{i,j - \frac{1}{2}} - F^a_{i,j + \frac{1}{2}} + F^d_{i,j - \frac{1}{2}} - F^d_{i,j + \frac{1}{2}}\Big) + S_{i,j}.
\end{array}
\f]
Here \f$F^a\f$ and \f$F^d\f$ represent the advective and diffusive fluxes across the cell boundaries. Since we assume ions and neutral particles to be stationary, these terms are nonzero only for electrons. For heavy particles, only the source term \f$S_{ij}\f$ remains.

The advective flux, \f$F^a\f$ uses an upwind scheme with flux limiting and is defined as follows:
\f[
 \begin{array}{ll}
F^a_{i + \frac{1}{2},j} = & E^+_{r; ~ i + \frac{1}{2},j} \Big[ n_{i,j} + \psi(P_{i,j})(n_{i+1,j} - n_{i,j}) \Big] \\
& E^-_{r; ~ i + \frac{1}{2},j} \Big[ n_{i + 1,j} + \psi(\frac{1}{P_{i+1,j}})(n_{i,j} - n_{i+1,j}) \Big]
\end{array}
\f]

\f[
 \begin{array}{ll}
   F^a_{i,j + \frac{1}{2}} = & E^+_{z; ~ i,j + \frac{1}{2}} \Big[ n_{i,j} + \psi(Q_{i,j})(n_{i,j+1} - n_{i,j}) \Big] \\
                             & E^-_{z; ~ i,j + \frac{1}{2}} \Big[ n_{i,j + 1} + \psi(\frac{1}{Q_{i,j+1}})(n_{i,j} - n_{i,j+1}) \Big],
\end{array}
\f]
where \f$E^+ = max(-E,0)\f$ and \f$E^- = min(-E,0)\f$ are used to distinguish the upwind directions for the components of the electric field, \f$E_r\f$ and \f$E_z\f$, and we have
\f[
\begin{array}{lll}
 P_{i,j} & = & \frac{n_{i,j} - n_{i-1,j}}{n_{i+1,j} - n_{i,j}}\\
 Q_{i,j} & = & \frac{n_{i,j} - n_{i,j-1}}{n_{i,j+1} - n_{i,j}}.
\end{array}
\f]
\f$\psi\f$ is the Koren limiter function:
\f[
 \psi(x) = max(0, min(1, \frac{1}{3} + \frac{x}{6}, x)).
\f]
The diffusive flux \f$F^d\f$ is calculated using a second-order central differences scheme:
\f[
\begin{array}{lll}
 F^d_{i + \frac{1}{2},j} & = & \frac{D}{\Delta r}(n_{i,j} - n_{i+1,j})\\
 F^d_{i,j + \frac{1}{2}} & = & \frac{D}{\Delta z}(n_{i,j} - n_{i,j+1})
\end{array}
\f]
and the reaction term \f$S_{i,j}\f$ is computed as
\f[
 S_{i,j} = \sum_{A~\in~{reactions}} \Big[ k_A(|\mathbf{E}|_{i,j}) \prod_{s~\in~{Spec(A)}} n_{s; i,j} \Big]
\f]
where \f$k_A\f$ denotes the field-dependent reaction rate coefficient of
reaction \f$A\f$, and \f$Spec(A)\f$ the set of species that appear as an input
for reaction \f$A\f$.

\subsubsection sect_dpe Discretization of the Poisson equation

We compute the net charge \f$q_{i,j}\f$ in a cell center by adding up the contributions from the individual charged species:
\f[
 q_{i,j} = \sum_{s~\in~{species}} n_{s; i,j} q_s.
\f]
With this net charge, the electric potential \f$\phi\f$ can be computed in the cell centers through a second-order central approximation of the dimensionless Poisson equation:
\f[
 q_{i,j} = \frac{\phi_{i+1,j} - 2 \phi_{i,j} + \phi_{i-1,j}}{\Delta r^2} + \frac{\phi_{i+1,j} - \phi_{i-1,j}}{2r_{i,j} \Delta r} + \frac{\phi_{i,j+1} - 2 \phi_{i,j} + \phi_{i,j-1}}{\Delta z^2}.
\f]
From the potential we can compute the components of the electric field from \f$\mathbf{E} = - \nabla \phi\f$ in the cell boundaries:
\f[
\begin{array}{lll}
 E_{r; ~ i + \frac{1}{2},j} & = & \frac{\phi_{i,j} - \phi_{i+1,j}}{\Delta r}\\
 E_{z; ~ i,j + \frac{1}{2}} & = & \frac{\phi_{i,j} - \phi_{i,j+1}}{\Delta r}.
\end{array}
\f]
The electric field strength is determined at the cell center, so we have to compute the field components in the center by averaging the values on the boundaries after which we can compute the field strength:
\f[
 |\mathbf{E}|_{i,j} = \sqrt{\left(\frac{E_{r;i - \frac{1}{2},j} + E_{r;i + \frac{1}{2},j}}{2}\right)^2 + \left(\frac{E_{z;i,j - \frac{1}{2}} + E_{z;i,j + \frac{1}{2}}}{2}\right)^2}.
\f]

\subsection sect_timestepping Time stepping

The code uses the explicit trapezoidal rule, a second order Runge-Kutta method,
for the temporal discretization with time step \f$\Delta t\f$.
Given some time step \f$t_i = i \Delta t\f$, density distributions
\f$\mathbf{n}_i(r,z) = \mathbf{n}(r,z,t_i)\f$ and electric field
\f$\mathbf{E}_i(r,z) = \mathbf{E}(r,z,t_i)\f$, the densities and field at the next
time step, \f$t_{i+1}\f$ are calculated by first computing an intermediate result
for the densities:
\f[
 \overline{\mathbf{n}}_{i+1} = \mathbf{n}_i + \Delta t F(\mathbf{n}_i, \mathbf{E}_i).
\f]
Using these intermediate densities, the potential can be computed by solving the
Poisson equation, after which we obtain the intermediate electric field
\f$\overline{\mathbf{E}}_{i+1}\f$.
With this, we compute the final values of the densities at \f$t_{i+1}\f$:
\f[
 \mathbf{n}_{i+1} = \mathbf{n}_i + \frac{\Delta t}{2} F(\mathbf{n}_i, \mathbf{E}_i) + \frac{\Delta t}{2} F(\overline{\mathbf{n}}_{i+1}, \overline{\mathbf{E}}_{i+1}).
\f]
Finally, we again compute the potential and electric field, now using the final
values of the densities.

The size of the time step \f$\Delta t\f$ is determined by using a Courant-Friederichs-Levy (CFL) restriction for stability of the advection part of the equations:
\f[
 \texttt{max} E_r \frac{\Delta t}{\Delta r} + \texttt{max} E_z \frac{\Delta_t}{\Delta_z} < \nu_a.
\f]
There are additional restrictions from other diffusion and reaction parts of the
equations, but they are dominated by the CFL criterior for the advection part, see
[Montijn et al.][1].
The value of \f$\nu_a\f$ is typically set to 0.25, which is well below the maximum
required for stability.
We refer the interested reader to the [Monograph by Hundsdorfer and Verwer][8].

\section sect_refinement Overview of refinement strategies and criteria

\subsection sect_overview Overview

The \c ARCoS simulation code contains functions for adaptive grid refinement (also known as adaptive mesh refinement or AMR). Since streamers span different length scales, there is a need to simulate relatively large physical domains while still having high spatial resolution in areas such as the streamer head. To ensure that such large domains can be simulated without giving up resolution and accuracy, the numerical grid is refined adaptively at each time step. The equations are solved on a coarse grid, after which the solution is analyzed using refinement criteria to determine the areas where refinement is needed. The equations are then solved on the refined subgrids after which the process is iterated. Grid generation and grid refinement are performed separately for the density equations and for the Poisson equation.

There are three main refinement criteria. The first two concern refinement of the density grids: Refinement based on the absolute value of \f$\mathbf{E}\f$ and refinement based on the curvature of densities (both charge density and particle density). The grids used by the [FISHPACK][9] solver use their own refinement scheme where the decision to refine is made if the difference between the solution on a grid and the solution on a finer grid exceeds a threshold. The [FISHPACK][9] solver is used for both the Poisson equation that determines the electric potential of the system and the Helmholtz equations for the photoionization reactions.

\subsection sect_size_refinement Size of the refined areas of the density grids

All CDR (Convection-Diffusion-Reaction, CDR is the shorthand term for the density part of the code) refinement criteria are on a per-point basis, which means that the question whether to refine or not is initially answered for every grid cell. This is inconvenient for several reasons, primarily due to the computational cost of such a scheme. The regions containing the streamer head will almost always need to be refined, it is not necessary to evaluate this point-by-point in these regions.

To ease this problem, a minimal refinement area is defined by two parameters:\n
\c cdr\_brick\_dr and \c cdr\_brick\_dz, see e.g., file default.cfg . The refinement module divides the grid it receives (this can be the coarsest grid covering the entire domain or a refined grid covering only part of the domain, the code and grid structure are recursive) into "bricks" of these dimensions and searches each brick for cells that match the refinement criteria. Once such a cell is found, the entire brick containing that cell is refined.

For the [FISHPACK][9] module, a different approach is used. The refinement routine scans its input grid, starting at the top (\f$z = z_{min}\f$), going down per "line" (a set of cells with equal \f$z\f$d-coordinate). Once it finds a line with points that meet the refinement criterion it searches for the first line that does not contain any points that meet the criterion. It then refines the smallest rectangular area that contains all the points that meet the criterion. This process is repeated until the bottom (\f$z = z_{max}\f$) of the grid is reached.

<img src="http://md-wiki.project.cwi.nl/images/refinement_cdr.png" width="360px" />
\par Figure 2. The nested structure of refined density grids

In <a name="refinement_cdr">Figure 2.</a>,
the black squares represent grid cells at the coarsest level (level 0),
the dark gray cells are the first refined sublevel (level 1). Two rectangular
grids are included at this level, their shared border is indicated by the
red line.
The light gray cells show grids at a further refined level (level 2).

<img src="http://md-wiki.project.cwi.nl/images/refinement_poisson.png" width="360px" />
\par Figure 3. The nested structure of refined Poisson grids

The black grid, as shown in <a name="refinement_poisson">Figure 3.</a> is the coarsest level, the dark gray cells are the first
refined sublevel, the light gray cells show grids at a further refined level.
Each grid has at most one subgrid.

The tree of grids for the density equations may contain refined grids that are adjacent to each other. A schematic showing the nested structure of refined density grids is shown in [Figure 2.](#refinement_cdr). The red line in this figure indicates the shared border between two subgrids. For the Poisson-grids, such a structure is not possible and a grid can have at most one refined child-grid as depicted in [Figure 3.](#refinement_poisson).

\subsection sect_criterion The |E| criterion

The electric field criterion is the most simple of the three refinement criteria. It is an empirical criterion that is not directly motivated by the underlying numerics. A cell with coordinates \f$(r,z)\f$ qualifies for refinement if:
\f[
|\mathbf{E}(r,z)| > E_c
\f]
where \f$E_c\f$ is the threshold electric field strength for refinement. \f$E_c\f$ is a user-determined parameter that is provided in the input file for a run. Since this criterion is independent of the grid level or the cell size, once a cell meets the criterion at the coarsest level, it will also do so at every refined level. Because of this property, the user can limit the refinement depth that is reached through this criterion with the \c ref\_level\_eabs input parameter, see e.g., file default.cfg. Setting \c ref\_level\_eabs to 1, for example, restricts the refinement from the coarsest level to the first refined level due to the \f$|\mathbf{E}|\f$ criterion.

The \f$|\mathbf{E}|\f$ criterion is inflexible in the sense that it requires the user to have advance knowledge of what the field strengths will be. A possible alternative would be to replace the fixed threshold value \f$E_c\f$ by a dimensionless fraction \f$c\f$ and refine if
\f[
|\mathbf{E}(r,z)| > c E_{max} 
\f]
with \f$E_{max}\f$ the maximum electric field strength in the computational domain. Since the electric field criterion is mostly empirical, picking the right value for the refinement threshold may be a trial-and-error process.

\subsection sect_curvature The curvature criteria

There are two criteria that use the curvature of density functions in order to determine which areas to refine. If the curvature is large compared to the cell size, the numerics may become unreliable and it is desirable to work with a finer grid. For a density function \f$u(r,z)\f$ and a cell size \f$\triangle r \times \triangle z\f$ the curvature function \f$C_u(r,z)\f$ is a discretization of the second derivative of \f$u\f$ in cylindrical coordinates \f$(r,z)\f$:
\f[
 \begin{array}{ll}
   C_u(r,z) = & \frac{1}{r + \frac{\triangle r}{2}}\Big[(r + \triangle r)\big(u(r + \triangle r, z) - u(r,z)\big) - r\big(u(r,z) - u(r - \triangle r,z)\big)\Big] + \\
              & \big[u(r, z + \triangle z) - 2 u(r,z) + u(r, z - \triangle z)\big].
 \end{array}
\f]
Rather than the absolute value of the curvature, the refinement module looks at the curvature relative to the global maximum, \f$Max(u)\f$. The final criterion then reads:\n
\n
\b Refine \f$(r,z)\f$ \b if \f$\frac{C_u(r,z)}{Max(u)} > C_t\f$\n
\n
with \f$C_t\f$ the threshold curvature. This refinement criterion is checked for two density functions \f$u\f$. The first is the (absolute) charge density function. Here an extra condition applies: the absolute value of the charge needs to exceed a certain threshold value (which is hard-coded) before a cell can qualify for refinement based on this criterion. Secondly, the curvature criterion is applied to the particle density functions. Since only mobile particles require a high spatial resolution, any immobile species are not considered in these criteria (which currently excludes all species other than electrons). The computational grids for these immobile species are simply the same as the grids used to solve the density equations for electrons.

\subsection sect_fishpack FISHPACK refinement

The [FISHPACK][9] module, for the Poisson equation and the photoionization equations, uses a different set of grids than the CDR module and with it a different refinement scheme. Initially, two grids are set up, one coarse and one fine grid (with the fine grid having twice the spatial resolution in each dimension, so 4 times the number of cells). The Poisson/Helmholtz equation is then solved on both grids and the solution of the coarse grid is interpolated onto the fine grid. A grid cell then qualifies for refinement if the absolute difference between the interpolated coarse solution and the fine solution (this difference is called the error) is more than some user-defined threshold. When refinement is needed, a new set of grids is determined using the strategy mentioned earlier and the process is repeated until either the desired accuracy is reached or the maximum number of allowed refinement levels is reached. Since the [FISHPACK][9] module was originally only used to solve the Poisson equation for the electrostatic problem and the value of the electric field 
is defined on the edge of a cell, a cell that does not meet the error-criterion still qualifies for refinement if its neighbor does meet the error-criterion.

One limitation to this scheme is the limited number of grid cells that the [FISHPACK][9] routine can handle. Since [FISHPACK][9] applies a cyclic reduction scheme, the round-off error increases with the number of grid cells. This places a limit on the size of grids that [FISHPACK][9] can solve. Once the refinement module wants to create a grid that is larger than the so-called [FISHPACK][9] limit, the refinement attempt is rejected and the code relaxes the error threshold by a factor of 2 and again determines the area to refine, using the new threshold.

To solve the photoionization problem, two Helmholtz equations need to be solved (For details on the implementation of photoionization, the reader is referred to [PhD Thesis Wormeester(to appear in 2013)][6], Chapter 3, section Photoionization 
and references therein). Each of the so-called "photo-terms" has its own characteristic absorption length, which depends on the gas density and oxygen ratio. The term with the short absorption length is often dominated by impact ionization in the head of the streamer, while the term with the long absorption length is the main contributor of electrons in front of the streamer head that are required for a positive streamer to propagate.

The default behavior of the \c ARCoS code is to treat these two photoionization terms in the same manner as the Poisson problem when it comes to refinement: all user-definable parameters were equal. Since the term with the short absorption length gives rise to a solution that benefits strongly from high spatial resolution (due to the steep gradients) it will easily trigger the refinement criterion. However, it is this term that is dominated by impact ionization, see [Luque et al.][2],
which reduces the relevance of accurate computation of this term. The user can therefore specify the refinement criteria for each of the two photoionization terms separately, providing the user with the means to allow the important, long absorption length term to benefit from high spatial resolution, while reducing the computational cost incurred by the less important term. However, in tests it was found that tuning the refinement criteria for the photoionization terms has very little effect on 
computational cost or results.

\subsection sect_Conclusions Conclusion
The adaptive refinement scheme of \c ARCoS allows for the simulation of large domains while maintaining high spatial resolution in regions that require this. A number of refinement parameters influence both the computational performance and the accuracy of the results, which means that the user has to monitor the results carefully. Since the refinement criteria were setup by
[Montijn et al.][1].
and
[Luque et al.][2].
for simulations of air and pure nitrogen, application of the code to other gases may require changes to the values of the various thresholds used in the refinement criteria. An example is high-purity oxygen, with a small nitrogen admixture. In such a gas, ionizing photons will have a very short characteristic absorption length and the calculation of the photoionization terms should be done with high accuracy close to the photon source, primarily the streamer head. However, the limitation of the [FISHPACK][9] refinement method does not permit several smaller, adjacent 
refined sub-grids, which makes it difficult to properly focus on the streamer head without including too much of the channel.

\section sect_software ARCoS software

\subsection sect_ARCoS_overview Basic overview and functionality

The \c ARCoS simulation software was originally developed by Alejendro Luque as a more flexible version of the adaptive refinement code developed by Carolynne Montijn as described in [PhD Thesis Montijn][7].
The original code by Montijn has been written in Fortran90, while \c ARCoS has been written in C.
The original [FISHPACK][9] package used for solving the Poisson and Helmholtz equations is written
in Fortran77 and was developed by Adams, Swarztrauber and Sweet. The \c ARCoS code is now compiled
with [FISH90][10], a modernization of the original [FISHPACK][9], employing Fortran90 to slightly
simplify and standardize the interface to some of the routines.

\c ARCoS solves the fluid equations for streamers, described in section \ref sect_fluid_eq,
on nested Cartesian grids using an adaptive mesh refinement technique.
\c ARCoS allows for the simulation of both positive and negative streamers in the
electrode configurations plate-plate and needle-plate.
The needle-plate electrode geometry is included using a charge simulation method
[Luque et al.][2].
This method replaces the electrode needle by a single point charge, with the
location and the size of the charge being updated at every time step to ensure the
potential at needle tip remains fixed at the predetermined value.
The limitation of this method is that the potential on the rest of the surface
of the simulated needle will not be accurate.
Consequently, the [continuity](#continuity) equation is only solved on
a smaller grid, not containing the simulated needle.

The effect of this is that \c ARCoS is not well suited for the study of the
inception of streamers, as the area around the tip of the needle is not
accurately modelled.
However, since inception is often affected by the behavior
of individual particles, the use of a particle code such as described in
[Teunissen et al.][5] and [Li et al.][4].
is recommended for studying streamer inception.
The purpose of the \c ARCoS code is to study streamer propagation in the phase
after the streamer has formed.
Studies performed by [Luque et al.][2] show that the dynamics of
streamers in later stages hardly depends on initial conditions.

\c ARCoS allows the user full control over the numerical parameters of the
simulation: grid size, refinement criteria and CFL numbers can be set by the user.
The kinetic model, i.e., the list of particle species, their reactions and
initial densities as well as the diffusion and mobility coefficients can be
specified via a series of input files, allowing the user to fine-tune the properties
of the gas in which the streamer is simulated, see configuration file input/kinetic_example.cfg.

The \c ARCoS code can be downloaded from the website
\c http://md-wiki.project.cwi.nl/

\subsection sect_IO Handling the software, input and output

\subsubsection sect_sim Starting a simulation

Two input parameter files governs all details of the simulation: 

\li \b Physical \b parameters such as voltage, electrode configuration, size of the gap, etc.
\li \b Numerical \b parameters such as grid size, refinement criteria, etc.
\li \b Practical \b details like the directory name, where the output files should be stored and the interval at which output should be generated.

[libconfig](http://www.hyperrealm.com/libconfig/), a free library for processing
structured configuration files, is used for reading, manipulating and writing these files.
The first file, stored as input/default.cfg, must contain the default values for the global variables.
This file is a part of the streamer package distribution.
The second file, say input/user_init.cfg, an example is given by input/example_user_init.cfg,
has a structure analogously to input/default.cfg, and should contain the parameters which differ
from the default values.
The program delivers a configuration file, say input/example_user_continue.cfg with the updated parameters
from input/user_init.cfg completed with the default values of input/default.cfg.

Since the execution time of a single run will take on the order of several days, it is recommended to
split the time period into smaller pieces, and restart the execution several times from the point
where the previous run stopped.
The easiest way to restart the simulation is
 \li  to copy the file input/user_continue.cfg into input/default.cfg, to be sure that equal values for the parameters are used, 
 \li to edit the file input/user_init.cfg, and change the \b \c t_end value, the \b \c restart value and the name of the \b \c load_file. See the end of of this section. 

The use of the configuration files construction has the following advantages:
 \li recompilation of the code is not necessary in case of a restart
 \li the user always has a clear overview of the parameter values used
 \li results of different or continuing runs can be stored in different output directories,
     as listed in the configuration file.
 \li besides the parameter value also comment coupled to a parameter can be changed in the configuration files written by the user. The length of the comment must be restricted by 100 characters
 \li the order of the parameters in the configuration file is free
 \li the user has the possibility to control the simulation, many parameter values can be changed.  

It is easy to resume a simulation by using a set of output data as initial conditions.
One has to adapt the parameter file with modified start and end times for the simulation.
To start a \c ARCoS simulation use the following command from
the directory containing the executables:
\n
\code
 ./arcos > out.example 2> err.example
\endcode
\n

The \c ARCoS program starts and it will print out the parameter values used:
 \li in \c out.example
 \li in input/user_continue.cfg.

The program will print some extra information to file \c out.example, e.g., the step size and 
when a new set of output data has been written, and to which set of file names.
Warnings and errors will be printed in file \c err.example. The program can terminate in three different ways:

 \li The preset end-time is reached.
 \li The program is terminated by the user.
 \li The time-step (as determined by the Courant criterion for stability, more details to come) has dropped below a preset threshold. This usually points to some form of instability.

In case the simulation runs on a PC or desktop machine, a convenient approach is to set a very
large value for the end time and, rather than having the program determine when to terminate,
keeping track of the progress of the streamer by checking the output files and manually terminating
the program when the desired output is reached, e.g., the streamer has reached the electrode, or,
 it has started to branch.
In other cases, it may be necessary for the program to be able to run for a fixed amount of time.
For example when it needs to exit gracefully, which is required for profiling software to work.
Also in case of batch jobbing with a limited CPU wall clock time, like on most supercomputers,
the end time must be chosen corresponding to the wall clock time.

Data files with periodical data controlled by \c $output_dt, stored in directory \c $output\_dir,
have names using the format \c variable.C123abc.tsv:
 \tparam <variable> is a particle species or electric field, e.g. \c charge, \c d_dummyminus, \c d_electrons \n
 \tparam <123> is the sequence number of that particular output dataset and \n
 \tparam <abc> specified the subgrid the output belongs to.

These files can be used to plot the solution or to restart the simulation.  More details on this in the section on output files. To resume the example simulation from output set \c 123, call:
\n
\code
./arcos example C123
\endcode
\n

\subsubsection sect_input The input files


File and directory handling:
 \li \b \c kin\_input - The filename of the [libconfig](http://www.hyperrealm.com/libconfig/)
			input file containing the species, seeds and reactions.
			By convention, these files have the extension \b .cfg.
 \li \b \c output\_dir - The directory where the output files will be stored.
			 This can be a relative path to the directory with the executable or an absolute
			path. When the execution starts, the directory $output_dir must be present and writable.
 \li \b \c output\_dt - The interval (in dimensionless units) with which an output dataset is to be saved.
			Lower values mean more frequent output, which gives finer grained time-dependent 
			data at the cost of more disk space.

Physical parameters:
 \li \b \c L\_r - Radius of the physical domain in dimensionless units.
 \li \b \c L\_z - Length of the physical domain. Does not include the needle.
 \li \b \c has\_photoionization - Whether to enable the photoionization module.\n 
                                  So, if \b \c has\_photoionization = 1, photoionization is present,\n
                                  if \b \c has\_photoionization = 0, execution without photoionization.\n
 \li \b \c photoionization\_file - Filename of the file containing the photoionization parameters.\n These parameters are further explained in appendix~ref{app:photoionization_parameters}.

\todo: WRITE APPENDIX

 \li \b \c E0\_z, \b \c E0\_y, \b \c E0\_x - Components of the external electric field in dimensionless units. Since the electrodes are located at the top and bottom of the domain \f$(z = L_z\f$ and \f$z = 0)\f$, only \c E0\_z should be nonzero. Positive values of \c E0\_z will generate an electric field in the \f$^z\f$ direction, with the top electrode (at \f$z = L_z\f$) having a negative charge, generating negative streamers and vice verse.
 \li \b \c pois\_inhom - Whether to use a needle-plane geometry or not.\n
                         If \b \c pois\_inhom = 1  then we have \b \c needle-plane case, \n
                         otherwise if \b \c  pois\_inhom = 0 gives the plane-plane case. \n
                         See below for additional remarks regarding the needle-plane geometry.
 \li \b \c needle\_length and \b \c needle\_radius - The length and radius of the needle in dimensionless units. Only applies when \b \c pois\_inhom = 1.
 \li \b \c end\_t - Time at which the simulation will stop, in dimensionless units.
 \li \b \c max\_ntheta - Number of azimuthal grid cells. Default \b \c max\_ntheta = 1. \n
                         \b \c max\_ntheta  > 1 will activate the full 3D simulation (without cylindrical symmetry) using a pseudo-spectral method described in more detail in [Luque et al.(2011)][3].

Numerical parameters:
 \li \b \c gridpoints\_r, \b \c gridpoints\_z - Number of gridpoints in \f$r\f$ and \f$z\f$ direction for the density equations at the coarsest level.
 \li \b \c cdr\_max\_level - Maximum number of refinement levels of the grid for the fluid equations. \n
       \b \c cdr\_max\_level = 0 means no refinement.
 \li \b \c cdr\_bnd\_bottom, \b \c cdr\_bnd\_top, \b \c cdr\_bnd\_right - Boundary condition for the density equations at the bottom \f$(z = 0)\f$, top \f$(z = L_z)\f$ and right \f$(r = L_r)\f$ of the domain. \n
 \b \c cdr\_bnd\_xxx = 1 gives homogeneous Neumann boundary conditions, \n
 \b \c cdr\_bnd\_xxx = -1 gives homogeneous Dirichlet boundary conditions.
 \li \b \c pois\_max\_level - Maximum number of refinement levels of the grid for the Poisson equation.
\b \c pois\_max\_level = 0 means no refinement.
 \li \b \c pois\_bnd\_bottom, \b \c pois\_bnd\_top, \b \c pois\_bnd\_right - Boundary condition for the Poisson equation at the bottom \f$(z = 0)\f$, top \f$(z = L_z)\f$ and right \f$(r = L_r)\f$ of the domain. \n
\b \c pois\_bnd\_xxx = 1 gives homogeneous Neumann boundary conditions, \n
\b \c pois\_bnd\_xxx = -1 gives homogeneous Dirichlet.

\subsubsection sect_parameter The parameter file

In the parameter file input/default.cfg variables and their values are assigned with the following syntax
\verbatim
{
    type = "type";       /* type can be a "string", a "double", an "int" */
    name = "name";       /* variable name as defined in file include/parameters.h */
    comment = "comment"; /* description of the variable, maximum of 100 characters */
    value = "value";     /* value of type string, double or integer, related to "type" */       
},

 # <-- this starts a comment-line, which is ignored.

 # White lines are also ignored.
 # String-values should be between quotes, numeric values should not
 pi=3.14
 # (Note: The above is an approximation)
 # Scientific notation can be used:
 pi_times_thousand=3.14E3
\endverbatim

Non-existent or misspelled parameters are ignored.
Parameters are not required, every parameter has a default value, which can be found
in function \b \c init_parameters of src/cstream.c. 

The following overview lists the parameters generally meant for testing purposes and changing them is not required for streamer simulations. Therefore they are best left at their default value.

File and directory handling:
 \li \b \c cdr\_output\_margin - Number of margin-cells to be added to the output of the density grids. 2 layers of ghost cells are added on the edge of the computational domain for the purpose of enforcing boundary conditions. With this parameter, these can be included in the output. Value must be smaller than or equal to 2. Default \b \c cdr\_output\_margin = 0.
 \li \b \c pois\_output - Output the grids used in solving the Poisson equation? The grids used for the Poisson equation are different than those used for the density equations as detailed in section \ref sect_refinement. \b \c pois\_output = 0. Note that the absolute electric field is already saved as a density grid, so it is not required to set this parameter to 1 to obtain this data.
 \li \b \c pois\_output\_margin - Number of margin-cells to be added to the output of the Poisson grids. Default  \b \c pois\_output\_margin = 0.  Note : \b \c pois\_output\_margin <= 2.

Physical parameters:
 \li \b \c start\_t - Initial time. This parameter is used when resuming simulations and is automatically updated in that case.

Numerical parameters: \n

The adaptive mesh refinement criteria and the parameters related to them are discussed in more detail in section \ref sect_refinement.

 \li \b \c ref\_threshold\_eabs - Refine grid if \f$|\mathbf{E}|\f$ exceeds this value.
 \li \b \c ref\_level\_eabs - Maximum number of refinement levels of the grid for the fluid equation due to the \f$|\mathbf{E}|\f$ criterion.
 \li \b \c ref\_threshold\_charge - Refinement threshold for the curvature of the charge density.
 \li \b \c ref\_threshold\_dens - Refinement threshold for the curvature of the particle densities.
 \li \b \c cdr\_brick\_r, \b \c cdr\_brick\_z - Size \f$(r\f$ and \f$z)\f$ of the minimal refinement area. See Section \ref sect_size_refinement.
 \li \b \c pois\_max\_error - An area of the Poisson grid is further refined if the relative error between two consecutive refined levels exceeds this value.
 \li \b \c nu\_a - Courant number based on advection to determine the time step. Note, that \b \c nu\_a < 1 to satisfy CFL stability. In streamer simulations, this time step restriction will dominate over other parameters (\b \c nu\_d (diffusion) and \b \c nu\_rt (reaction)). More details on the time stepping can be found in \ref sect_timestepping.

\subsubsection sect_needle Implementation of the needle-plane configuration

As mentioned in section \ref sect_overview , the needle-plane electrode geometry is
implemented using a charge simulation technique, which means that the entire needle is
represented by a single point charge located on the axis.
The position and strength of this charge is updated each time step to ensure that the
potential at the point that would be the tip of the needle remains fixed.
A schematic depiction of this setup can be seen in [Figure 1.](#setup).
This provides a reasonable approximation of the potential in the area below
the needle (the needle is always located at the top of the domain), but the potential
will be wrong in the areas to the sides of the needle.
Consequently, the density equations are only solved from the tip of the needle and downwards.
This means that the computational domain for the density equations is smaller than
that for the Poisson equation.

In the specification of the external electric field, \c E0\_z is interpreted as the electric field between the two planar electrodes, far away from the needle. The applied potential is computed as follows:
\f[
  V = E_{0,z} * (L_z + L_{needle})
\f]

\subsubsection sect_output Output

The code generates a large amount of output files at every \c \b output\_dt units
of simulated time.
These output files combined are sufficient to resume the simulation
from that time step and can also be used for data analysis.
Output files are named \c variable.C123abc.tsv, where 'tsv' stands for 'Tab-Separated Values'.
\c variable is a particle species or electric field.
Examples include \c electrons, \c n2plus and \c eabs
(the latter being \f$|\mathbf{E}|\f$). \c 123 is the sequence number or
output time step (which is unrelated to the time step of the numerical scheme) of that
output set, it starts at \c 000 for the first set.
This is the initial condition of the system before the simulation starts at \f$t = 0\f$.
The alphabetic part of the filename after the output time step, \c abc denotes the
subgrid contained in that file and this extension is defined recursively.
The coarsest grid covering the entire domain has no such alphabetical extension
(example: \c electrons.C123.tsv).
At every next level, the subgrids at that level are named \c a, \c b, \c c, ...
and this letter is appended to the alphabetical extension.
The file \c electrons.C123ba.tsv contains electron densities at the
123-th time step (so at \f$t = 123 * output\_dt\f$ ) for the first
subgrid of the second subgrid of the main grid.
Beware that since the grid refinement is adaptive, each time step has different subgrids.

<img src="http://md-wiki.project.cwi.nl/images/grids.png" width="360px" />
\par Figure 4. Schematic depiction of the naming convention for output files. Each next letter corresponds to a new refined level.

Data, as shown in <a name="output_grids">Figure 4.</a>, is stored as plain text,
with each line containing a single number.
Data is ordered in columns (with fixed \f$r\f$ coordinate).
So to read the data, use the following pseudo-code:
\verbatim
for (i = 0, i < rmax*zmax; i++)
{
  r = floor(i / zmax);
  z = i % zmax;
  data[r,z] = read_line_from_file();
}
\endverbatim

The dimensions of the grid are not contained in the data files.
Instead, for each subgrid and output-step, two additional files are created: 

1. _r.C123abc.tsv_\n
2. _z.C123abc.tsv_\n
\n

corresponding to subgrid _abc_ of output step _123_. The structure of these files is identical to that of the regular data files, but instead of particle densities or electric field strengths, these files contain the \f$r\f$ and \f$z\f$ coordinates of the center of the cell corresponding to that line-number. So to determine the coordinates of the \f$n^{th}\f$ line in a regular data file, simply read the \f$n^{th}\f$ line of the corresponding \b r and \b z files.

\subsubsection sect_structure File structure

<img src="http://md-wiki.project.cwi.nl/images/structure.png" width="500px" />
\par Figure 5. File structure of of ARCoS package.

The distribution  of the ARCoS package contains several directories:

\li \b \c FISH90 - the directory where FISH90 or FISHPACK should be downloaded

\li \b \c output - this directory may be empty. Its name must correspond to the value of \b \c output_dir in file input/default.cfg. Or, the value in input/user_init.cfg if present.
\li \b \c python - contains plot files. By means of \b \c plotvar pictures can be made of the output files. Not yet implemented. 
\li \b \c doc - directory with files for [Doxygen][12] to generate documentation from source code and this \b \c MANUAL. The source of this manual (in MANUAL.md) can be found in its directory \b \c markdown. See also the \b \c doxygen_config_file.
\li \b \c arcos_f90 - a part of the \b \c functions in file cdr.c have been replaced by a piece of \b \c FORTRAN90 code in order to accelerate the simulation. 
\li \b \c src - this directory contains the source files written in \b \c c.
\li \b \c include - this directory contains the include files
\li \b \c input - this directory contains all input files. Most of them are libconfig configuration files.


\subsubsection sect_source Source files

The \c ARCoS simulation software was mostly written in \c C and its source code is split up in several files, each dealing with a separate part of the program. Most source files have an associated header file (the source file \b \c example.c  has header file \b \c example.h ) containing the type-definitions and preprocessor macros. The function prototypes are aggregated in the header file \c proto.h. Below is a short summary of the important source files and the functionality that is contained within them.

 \li cdr.c - Functions for solving the convection-diffusion-reaction (CDR) equations, creation, manipulation and refinement of CDR grids and time stepping.
 \li configuration.c - Module for input/output of parameters. The code uses the library
[libconfig](http://www.hyperrealm.com/libconfig/)
 \li cstream.c - Contains some general initialization and termination functions.
 \li dft.c - Functions related with discrete fourier transformations.
 \li grid.c - Low-level functions for handling of grids, both CDR and Poisson grids.
 \li interpol2.c - Interpolation functions for the mapping of one grid to another (for example during refinement).
 \li main.c - Functions for reading input parameters, starting of the code and the main loop.
 \li mapper.c - Mapping of one grid tree onto another.
 \li misc.c - Miscellaneous utilities for allocating memory.
 \li photo.c - Photoionization functions.
 \li poisson.c - Functions for solving the Poisson equation, including manipulation of Poisson grids and calling the external [FISPACK][9] solver.
 \li reaction.c - Functions for computation of reactions between species as part of the density equations.
 \li react_table.c - Performs initialization of reaction coefficient tables as well as table lookups.
 \li rt.c - Functions for handling the loading of the input file containing the kinetic model (species, reactions, seeds).
 \li rz_array.c - Low-level functions for handling Fortran-compatible arrays.
 \li sprites.c - Routines for the sprites module.
