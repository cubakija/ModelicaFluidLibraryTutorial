package ShallowWellPumpingSystem
  model StopValve "Valve for (almost) incompressible fluids"
    extends Modelica.Fluid.Valves.BaseClasses.PartialValve(CvData = Modelica.Fluid.Types.CvTypes.Av, Av = 3.68008e-5, dp_nominal = 1e5, m_flow_nominal = 0.3);
    import Modelica.Fluid.Types.CvTypes;
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    import Modelica.Fluid.Utilities;
    constant SI.ReynoldsNumber Re_turbulent = 4000 "cf. straight pipe for fully open valve -- dp_turbulent increases for closing valve";
    parameter Boolean use_Re = system.use_eps_Re "= true, if turbulent region is defined by Re, otherwise by m_flow_small" annotation(
      Dialog(tab = "Advanced"),
      Evaluate = true);
    //SI.MassFlowRate m_flow_turbulent=if not use_Re then m_flow_small else
    //  max(m_flow_small,
    //      (Modelica.Constants.pi/8)*sqrt(max(relativeFlowCoefficient,0.001)*Av*4/pi)*(Medium.dynamicViscosity(state_a) + Medium.dynamicViscosity(state_b))*Re_turbulent);
    //SI.AbsolutePressure dp_turbulent_=if not use_Re then dp_small else
    //  max(dp_small, m_flow_turbulent^2/(max(relativeFlowCoefficient,0.001)^2*Av^2*(Medium.density(state_a) + Medium.density(state_b))/2));
    // substitute m_flow_turbulent into dp_turbulent
    SI.AbsolutePressure dp_turbulent = if not use_Re then dp_small else max(dp_small, (Medium.dynamicViscosity(state_a) + Medium.dynamicViscosity(state_b)) ^ 2 * pi / 8 * Re_turbulent ^ 2 / (max(relativeFlowCoefficient, 0.001) * Av * (Medium.density(state_a) + Medium.density(state_b))));
    Modelica.Blocks.Tables.CombiTable1D valveCharacteristicTable(smoothness = Modelica.Blocks.Types.Smoothness.MonotoneContinuousDerivative1, table = [0, 0; 0.052631579, 0.271493213; 0.157894737, 0.542986425; 0.210526316, 0.601809955; 0.263157895, 0.633484163; 0.368421053, 0.791855204; 0.473684211, 0.882352941; 0.631578947, 0.972850679; 0.789473684, 0.995475113; 1, 1]) annotation(
      Placement(visible = true, transformation(origin = {106, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  protected
    Real relativeFlowCoefficient;
  initial equation
    if CvData == CvTypes.OpPoint then
      m_flow_nominal = valveCharacteristic(opening_nominal) * Av * sqrt(rho_nominal) * Utilities.regRoot(dp_nominal, dp_small) "Determination of Av by the operating point";
    end if;
  equation
    connect(valveCharacteristicTable.u[1], opening_actual) annotation(
      Line(points = {{94, 20}, {70, 20}}, color = {0, 0, 127}));
// m_flow = valveCharacteristic(opening)*Av*sqrt(d)*sqrt(dp);
    relativeFlowCoefficient = valveCharacteristicTable.y[1];
    if checkValve then
      m_flow = homotopy(relativeFlowCoefficient * Av * sqrt(Medium.density(state_a)) * Utilities.regRoot2(dp, dp_turbulent, 1.0, 0.0, use_yd0 = true, yd0 = 0.0), relativeFlowCoefficient * m_flow_nominal * dp / dp_nominal);
/* In Modelica 3.1 (Disadvantage: Unnecessary event at dp=0, and instead of smooth=2)
  m_flow = valveCharacteristic(opening)*Av*sqrt(Medium.density(state_a))*
                (if dp>=0 then Utilities.regRoot(dp, dp_turbulent) else 0);
  */
    elseif not allowFlowReversal then
      m_flow = homotopy(relativeFlowCoefficient * Av * sqrt(Medium.density(state_a)) * Utilities.regRoot(dp, dp_turbulent), relativeFlowCoefficient * m_flow_nominal * dp / dp_nominal);
    else
      m_flow = homotopy(relativeFlowCoefficient * Av * Utilities.regRoot2(dp, dp_turbulent, Medium.density(state_a), Medium.density(state_b)), relativeFlowCoefficient * m_flow_nominal * dp / dp_nominal);
/* In Modelica 3.1 (Disadvantage: Unnecessary event at dp=0, and instead of smooth=2)
  m_flow = smooth(0, Utilities.regRoot(dp, dp_turbulent)*(if dp>=0 then sqrt(Medium.density(state_a)) else sqrt(Medium.density(state_b))));
  */
    end if;
    annotation(
      Documentation(info = "<html>
<p>
Valve model according to the IEC 534/ISA S.75 standards for valve sizing, incompressible fluids.</p>

<p>
The parameters of this model are explained in detail in
<a href=\"modelica://Modelica.Fluid.Valves.BaseClasses.PartialValve\">PartialValve</a>
(the base model for valves).
</p>

<p>
This model assumes that the fluid has a low compressibility, which is always the case for liquids.
It can also be used with gases, provided that the pressure drop is lower than 0.2 times the absolute pressure at the inlet, so that the fluid density does not change much inside the valve.</p>

<p>
If <code>checkValve</code> is false, the valve supports reverse flow, with a symmetric flow characteristic curve. Otherwise, reverse flow is stopped (check valve behaviour).
</p>

<p>
The treatment of parameters <strong>Kv</strong> and <strong>Cv</strong> is
explained in detail in the
<a href=\"modelica://Modelica.Fluid.UsersGuide.ComponentDefinition.ValveCharacteristics\">User's Guide</a>.
</p>

</html>", revisions = "<html>
<ul>
<li><em>2 Nov 2005</em>
  by <a href=\"mailto:francesco.casella@polimi.it\">Francesco Casella</a>:<br>
     Adapted from the ThermoPower library.</li>
</ul>
</html>"));
  end StopValve;

  model StopValveTest1
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.StopValve stopValve1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {0, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325 + 1e5) annotation(
      Placement(visible = true, transformation(origin = {-62, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {60, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10) annotation(
      Placement(visible = true, transformation(origin = {-30, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {58, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ramp1.y, stopValve1.opening) annotation(
      Line(points = {{-18, 32}, {0, 32}, {0, 4}, {0, 4}}, color = {0, 0, 127}));
    connect(stopValve1.port_b, boundary1.ports[1]) annotation(
      Line(points = {{10, -4}, {50, -4}, {50, -6}, {50, -6}}, color = {0, 127, 255}));
    connect(boundary.ports[1], stopValve1.port_a) annotation(
      Line(points = {{-52, -6}, {-10, -6}, {-10, -4}, {-10, -4}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"),
  __OpenModelica_commandLineOptions = "--matchingAlgorithm=PFPlusExt --indexReductionMethod=dynamicStateSelection -d=initialization,NLSanalyticJacobian  -d=aliasConflicts ");
  end StopValveTest1;

  model StopValveTest2
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.StopValve stopValve1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {32, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {72, -10}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 1, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {6, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {70, 28}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.OpenTank tank(redeclare package Medium = Medium, crossArea = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, height = 40, level_start = 35, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 1, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.015, height = 0.08)}, use_portsData = true) annotation(
      Placement(visible = true, transformation(origin = {-64, -16}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.015, height_ab = 2, length = 2) annotation(
      Placement(visible = true, transformation(origin = {-28, -26}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-28, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{-28, -16}, {-28, -10}}, color = {0, 127, 255}));
    connect(volume.ports[2], stopValve1.port_a) annotation(
      Line(points = {{-28, -10}, {22, -10}}, color = {0, 127, 255}));
    connect(tank.ports[1], pipe.port_a) annotation(
      Line(points = {{-64, -36}, {-60, -36}, {-60, -46}, {-28, -46}, {-28, -36}, {-28, -36}}, color = {0, 127, 255}));
    connect(ramp1.y, stopValve1.opening) annotation(
      Line(points = {{17, 28}, {32, 28}, {32, -2}}, color = {0, 0, 127}));
    connect(stopValve1.port_b, boundary1.ports[1]) annotation(
      Line(points = {{42, -10}, {62, -10}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 1200, Tolerance = 1e-06, Interval = 1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end StopValveTest2;

  function wellPumpHead
    extends Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.baseFlow;
  algorithm
    head := (-70000.0 * V_flow) + 47.0;
  end wellPumpHead;

  function wellPumpPower
    extends Modelica.Fluid.Machines.BaseClasses.PumpCharacteristics.basePower;
  algorithm
    consumption := 250000 * V_flow + 100;
  end wellPumpPower;

  model WellPump
    extends Modelica.Fluid.Machines.PrescribedPump(redeclare replaceable package Medium = Modelica.Media.Water.StandardWater, redeclare function flowCharacteristic = wellPumpHead, redeclare function powerCharacteristic = wellPumpPower, N_nominal = 3000, checkValve = true, use_N_in = true, m_flow_start = 0.4, use_powerCharacteristic = true, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.SteadyStateInitial);
  equation

    annotation(
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end WellPump;

  model WellPumpTest1
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.WellPump wellPump1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-38, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-82, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, nPorts = 1, use_p_in = true) annotation(
      Placement(visible = true, transformation(origin = {46, 4}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.025, length = 1) annotation(
      Placement(visible = true, transformation(origin = {2, 6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 400000, offset = 101325, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {44, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 3000) annotation(
      Placement(visible = true, transformation(origin = {-58, 44}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {86, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(const.y, wellPump1.N_in) annotation(
      Line(points = {{-47, 44}, {-38, 44}, {-38, 16}}, color = {0, 0, 127}));
    connect(ramp1.y, boundary1.p_in) annotation(
      Line(points = {{55, 44}, {70, 44}, {70, 12}, {58, 12}}, color = {0, 0, 127}));
    connect(pipe.port_b, boundary1.ports[1]) annotation(
      Line(points = {{12, 6}, {36, 6}, {36, 4}, {36, 4}}, color = {0, 127, 255}));
    connect(wellPump1.port_b, pipe.port_a) annotation(
      Line(points = {{-28, 6}, {-6, 6}, {-6, 6}, {-8, 6}}, color = {0, 127, 255}));
    connect(boundary.ports[1], wellPump1.port_a) annotation(
      Line(points = {{-72, 6}, {-48, 6}, {-48, 6}, {-48, 6}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 10, Tolerance = 1e-06, Interval = 0.02),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "ida", nls = "hybrid"));
  end WellPumpTest1;

  model WellPumpTest2
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.StopValve stopValve1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {90, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {130, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 200, height = 0.9998, offset = 0.0001, startTime = 2000) annotation(
      Placement(visible = true, transformation(origin = {62, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {128, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.OpenTank tank(redeclare package Medium = Medium, crossArea = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, height = 50, level_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.025, height = 0.08), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.015, height = 0.08)}, use_portsData = true) annotation(
      Placement(visible = true, transformation(origin = {-6, -12}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.015, height_ab = 2, length = 2) annotation(
      Placement(visible = true, transformation(origin = {30, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {30, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShallowWellPumpingSystem.WellPump wellPump1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-86, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe1(redeclare package Medium = Medium, diameter = 0.025, height_ab = 0, length = 1) annotation(
      Placement(visible = true, transformation(origin = {-52, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-124, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 3000) annotation(
      Placement(visible = true, transformation(origin = {-106, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{30, -12}, {30, -6}}, color = {0, 127, 255}));
    connect(tank.ports[2], pipe.port_a) annotation(
      Line(points = {{-6, -32}, {-2, -32}, {-2, -42}, {30, -42}, {30, -32}}, color = {0, 127, 255}));
    connect(const.y, wellPump1.N_in) annotation(
      Line(points = {{-94, 12}, {-86, 12}, {-86, -10}, {-86, -10}}, color = {0, 0, 127}));
    connect(ramp1.y, stopValve1.opening) annotation(
      Line(points = {{73, 32}, {90, 32}, {90, 2}}, color = {0, 0, 127}));
    connect(pipe1.port_b, tank.ports[1]) annotation(
      Line(points = {{-42, -20}, {-34, -20}, {-34, -42}, {-8, -42}, {-8, -32}, {-6, -32}}, color = {0, 127, 255}));
    connect(boundary.ports[1], wellPump1.port_a) annotation(
      Line(points = {{-114, -20}, {-96, -20}}, color = {0, 127, 255}));
    connect(wellPump1.port_b, pipe1.port_a) annotation(
      Line(points = {{-76, -20}, {-69, -20}, {-69, -20}, {-62, -20}, {-62, -20}, {-62, -20}}, color = {0, 127, 255}));
    connect(volume.ports[2], stopValve1.port_a) annotation(
      Line(points = {{30, -6}, {80, -6}}, color = {0, 127, 255}));
    connect(stopValve1.port_b, boundary1.ports[1]) annotation(
      Line(points = {{100, -6}, {120, -6}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 5000, Tolerance = 1e-06, Interval = 1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "ida", nls = "hybrid"),
      Diagram(coordinateSystem(extent = {{-160, -100}, {160, 100}})),
      Icon(coordinateSystem(extent = {{-160, -100}, {160, 100}})),
      __OpenModelica_commandLineOptions = "");
  end WellPumpTest2;

  model WellPumpTest3
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.StopValve stopValve1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {90, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {130, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 180, height = 1, offset = 0, startTime = 3000) annotation(
      Placement(visible = true, transformation(origin = {62, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {128, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Vessels.OpenTank tank(redeclare package Medium = Medium, crossArea = 0.01, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, height = 50, level_start = 0, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.025, height = 0.08), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.015, height = 0.08)}, use_portsData = true) annotation(
      Placement(visible = true, transformation(origin = {-6, -12}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.015, height_ab = 2, length = 2) annotation(
      Placement(visible = true, transformation(origin = {30, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {30, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShallowWellPumpingSystem.WellPump wellPump1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-86, -20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe1(redeclare package Medium = Medium, diameter = 0.025, height_ab = -0.22, length = 0.22) annotation(
      Placement(visible = true, transformation(origin = {-52, -30}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-160, -68}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 3000) annotation(
      Placement(visible = true, transformation(origin = {-106, 12}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe2(redeclare package Medium = Medium, diameter = 0.06, height_ab = 6.3, length = 6.3) annotation(
      Placement(visible = true, transformation(origin = {-136, -42}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume1(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start = system.p_start - 1000 * 6.3 * Modelica.Constants.g_n, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-136, -10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(boundary.ports[1], pipe2.port_a) annotation(
      Line(points = {{-150, -68}, {-136, -68}, {-136, -52}}, color = {0, 127, 255}));
    connect(pipe2.port_b, volume1.ports[1]) annotation(
      Line(points = {{-136, -32}, {-136, -20}}, color = {0, 127, 255}));
    connect(tank.ports[2], pipe.port_a) annotation(
      Line(points = {{-6, -32}, {-2, -32}, {-2, -52}, {30, -52}, {30, -32}}, color = {0, 127, 255}));
    connect(pipe1.port_b, tank.ports[1]) annotation(
      Line(points = {{-52, -40}, {-52, -52}, {-8, -52}, {-8, -32}, {-6, -32}}, color = {0, 127, 255}));
    connect(wellPump1.port_b, pipe1.port_a) annotation(
      Line(points = {{-76, -20}, {-52, -20}}, color = {0, 127, 255}));
    connect(volume1.ports[2], wellPump1.port_a) annotation(
      Line(points = {{-136, -20}, {-96, -20}}, color = {0, 127, 255}));
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{30, -12}, {30, -6}}, color = {0, 127, 255}));
    connect(const.y, wellPump1.N_in) annotation(
      Line(points = {{-94, 12}, {-86, 12}, {-86, -10}, {-86, -10}}, color = {0, 0, 127}));
    connect(ramp1.y, stopValve1.opening) annotation(
      Line(points = {{73, 32}, {90, 32}, {90, 2}}, color = {0, 0, 127}));
    connect(volume.ports[2], stopValve1.port_a) annotation(
      Line(points = {{30, -6}, {80, -6}}, color = {0, 127, 255}));
    connect(stopValve1.port_b, boundary1.ports[1]) annotation(
      Line(points = {{100, -6}, {120, -6}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 6000, Tolerance = 1e-06, Interval = 1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "ida", nls = "hybrid"),
      Diagram(coordinateSystem(extent = {{-200, -100}, {160, 100}})),
      Icon(coordinateSystem(extent = {{-200, -100}, {160, 100}})),
      __OpenModelica_commandLineOptions = "");
  end WellPumpTest3;

  model PressureTank "Simple tank with inlet/outlet ports"
    import Modelica.Constants.pi;
    import SI = Modelica.SIunits;
    import Modelica.Fluid.Types;
    import Modelica.Fluid;
    // Air
    replaceable package Air = Modelica.Media.Air.DryAirNasa;
    Air.BaseProperties air;
    SI.Volume V_air;
    SI.Mass M_air;
    SI.Energy U_air;
    // Tank properties
    SI.Height level(stateSelect = StateSelect.prefer, start = level_start_eps) "Level height of tank";
    SI.Volume V(stateSelect = StateSelect.never) "Actual tank volume";
    // Tank geometry
    parameter SI.Height height = 0.2 "Height of tank";
    parameter SI.Area crossArea = 0.1 "Area of tank";
    parameter SI.Volume V_total = height * crossArea;
    // Ambient
    parameter Medium.AbsolutePressure p_ambient = system.p_ambient "Tank surface pressure" annotation(
      Dialog(tab = "Assumptions", group = "Ambient"));
    parameter Medium.Temperature T_ambient = system.T_ambient "Tank surface Temperature" annotation(
      Dialog(tab = "Assumptions", group = "Ambient"));
    // Initialization
    parameter SI.Height level_start(min = 0) = 0.08 "Start value of tank level" annotation(
      Dialog(tab = "Initialization"));
    // Mass and energy balance, ports
    extends Modelica.Fluid.Vessels.BaseClasses.PartialLumpedVessel(final fluidVolume = V, final fluidLevel = level, final fluidLevel_max = height, final vesselArea = crossArea, heatTransfer(surfaceAreas = {crossArea + 2 * sqrt(crossArea * pi) * level}), final initialize_p = false, final p_start = p_ambient, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, portsData = {Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.025, height = 0.08), Modelica.Fluid.Vessels.BaseClasses.VesselPortsData(diameter = 0.015, height = 0.08)}, use_portsData = true);
    Modelica.Blocks.Interfaces.RealOutput tankPressure annotation(
      Placement(visible = true, transformation(origin = {60, 50}, extent = {{-10, -10}, {10, 10}}, rotation = 0), iconTransformation(origin = {-78, 110}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
  protected
    final parameter SI.Height level_start_eps = max(level_start, Modelica.Constants.eps);
  equation
// Total quantities
    V = crossArea * level "Volume of fluid";
    medium.p = air.p;
    tankPressure = air.p - p_ambient;
    V_air = V_total - V;
    M_air = air.d * V_air;
    U_air = air.u * M_air;
    der(M_air) = 0;
// Source termsEnergy balance
    if Medium.singleState or energyDynamics == Types.Dynamics.SteadyState then
      Wb_flow = 0 "Mechanical work is neglected, since also neglected in medium model (otherwise unphysical small temperature change, if tank level changes)";
      der(U_air) = 0;
    else
      Wb_flow = -p_ambient * der(V);
      der(U_air) = -Wb_flow;
    end if;
//Determine port properties
    for i in 1:nPorts loop
      vessel_ps_static[i] = max(0, level - portsData_height[i]) * system.g * medium.d + air.p;
    end for;
  initial equation
    if massDynamics == Types.Dynamics.FixedInitial then
      level = level_start_eps;
    elseif massDynamics == Types.Dynamics.SteadyStateInitial then
      der(level) = 0;
    end if;
    air.p = p_ambient;
    air.T = T_ambient;
    medium.p = p_ambient;
    annotation(
      defaultComponentName = "tank",
      Icon(coordinateSystem(initialScale = 0.2), graphics = {Rectangle(fillColor = {85, 170, 255}, fillPattern = FillPattern.VerticalCylinder, extent = {{-100, -100}, {100, 10}}), Text(extent = {{-95, 60}, {95, 40}}, textString = "level ="), Text(extent = {{-95, -24}, {95, -44}}, textString = "%level_start", fontName = "0"), Rectangle(origin = {-27, 60}, lineColor = {255, 170, 255}, fillColor = {255, 255, 255}, fillPattern = FillPattern.VerticalCylinder, extent = {{-73, 40}, {127, -50}}), Rectangle(origin = {-52, 64}, extent = {{-48, 36}, {152, -132}})}),
      Documentation(info = "<html>
<p>
Model of a tank that is open to the ambient at the fixed pressure
<code>p_ambient</code>.
</p>
<p>
The vector of connectors <strong>ports</strong> represents fluid ports at configurable heights, relative to the bottom of tank.
Fluid can flow either out of or in to each port.
</p>
The following assumptions are made:
<ul>
<li>The tank is filled with a single or multiple-substance medium having a density higher than the density of the ambient medium.</li>
<li>The fluid has uniform density, temperature and mass fractions</li>
<li>No liquid is leaving the tank through the open top; the simulation breaks with an assertion if the liquid level growths over the height.</li>
</ul>
<p>
The port pressures represent the pressures just after the outlet (or just before the inlet) in the attached pipe.
The hydraulic resistances <code>portsData.zeta_in</code> and <code>portsData.zeta_out</code> determine the dissipative pressure drop between tank and port depending on
the direction of mass flow. See <a href=\"modelica://Modelica.Fluid.Vessels.BaseClasses.VesselPortsData\">VesselPortsData</a> and <em>[Idelchik, Handbook of Hydraulic Resistance, 2004]</em>.
</p>
<p>
With the setting <code>use_portsData=false</code>, the port pressure represents the static head
at the height of the respective port.
The relationship between pressure drop and mass flow rate at the port must then be provided by connected components;
Heights of ports as well as kinetic and potential energy of fluid entering or leaving are not taken into account anymore.
</p>
</html>", revisions = "<html>
<ul>
<li><em>Dec. 12, 2008</em> by R&uuml;diger Franke: move port definitions
   to BaseClasses.PartialLumpedVessel; also use energy and mass balance from common base class</li>
<li><em>Dec. 8, 2008</em> by Michael Wetter (LBNL):<br>
Implemented trace substances.</li>
<li><em>Jan. 6, 2006</em> by Katja Poschlad, Manuel Remelhe (AST Uni Dortmund),
   Martin Otter (DLR):<br>
   Implementation based on former tank model.</li>
<li><em>Oct. 29, 2007</em> by Carsten Heinrich (ILK Dresden):<br>
Adapted to the new fluid library interfaces:
<ul> <li>FluidPorts_b is used instead of FluidPort_b (due to it is defined as an array of ports)</li>
    <li>Port name changed from port to ports</li></ul>Updated documentation.</li>
<li><em>Apr. 25, 2006</em> by Katrin Pr&ouml;l&szlig; (TUHH):<br>
Limitation to bottom ports only, added inlet and outlet loss factors.</li>
</ul>
</html>"),
      experiment(StartTime = 0, StopTime = 1, Tolerance = 1e-6, Interval = 0.002),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "dassl"));
  end PressureTank;

  model PressureTankTest1
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    PressureTank tank(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {2, 0}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary(redeclare package Medium = Medium, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {-44, -40}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.MassFlowSource_T boundary1(redeclare package Medium = Medium, m_flow = 0, nPorts = 1, use_m_flow_in = true) annotation(
      Placement(visible = true, transformation(origin = {50, -40}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 0.1, offset = 0, startTime = 0) annotation(
      Placement(visible = true, transformation(origin = {-84, -4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {68, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp2(duration = 10, height = -0.2, offset = 0, startTime = 120) annotation(
      Placement(visible = true, transformation(origin = {68, 0}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(ramp2.y, boundary1.m_flow_in) annotation(
      Line(points = {{80, 0}, {88, 0}, {88, -32}, {60, -32}, {60, -32}, {60, -32}}, color = {0, 0, 127}));
    connect(ramp1.y, boundary.m_flow_in) annotation(
      Line(points = {{-72, -4}, {-62, -4}, {-62, -32}, {-54, -32}, {-54, -32}}, color = {0, 0, 127}));
    connect(tank.ports[2], boundary1.ports[1]) annotation(
      Line(points = {{2, -20}, {6, -20}, {6, -40}, {40, -40}}, color = {0, 127, 255}, thickness = 0.5));
    connect(boundary.ports[1], tank.ports[1]) annotation(
      Line(points = {{-34, -40}, {-2, -40}, {-2, -20}, {2, -20}}, color = {0, 127, 255}, thickness = 0.5));
    annotation(
      experiment(StartTime = 0, StopTime = 250, Tolerance = 1e-6, Interval = 0.5),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", nls = "hybrid", outputFormat = "mat", s = "ida"));
  end PressureTankTest1;

  model PressureTankTest2
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.StopValve stopValve1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {90, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {130, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Ramp ramp1(duration = 10, height = 1, offset = 0, startTime = 35) annotation(
      Placement(visible = true, transformation(origin = {62, 32}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-160, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShallowWellPumpingSystem.PressureTank tank(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-46, -10}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.015, height_ab = 2, length = 2) annotation(
      Placement(visible = true, transformation(origin = {42, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {42, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShallowWellPumpingSystem.WellPump wellPump1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-106, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe1(redeclare package Medium = Medium, diameter = 0.025, height_ab = -0.22, length = 0.22) annotation(
      Placement(visible = true, transformation(origin = {-78, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-160, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Constant const(k = 3000) annotation(
      Placement(visible = true, transformation(origin = {-126, 56}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe2(redeclare package Medium = Medium, diameter = 0.06, height_ab = 6.3, length = 6.3) annotation(
      Placement(visible = true, transformation(origin = {-136, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume1(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start = system.p_start - 1000 * 6.3 * Modelica.Constants.g_n, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-136, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{42, -12}, {42, -6}}, color = {0, 127, 255}));
    connect(tank.ports[2], pipe.port_a) annotation(
      Line(points = {{-46, -30}, {-40, -30}, {-40, -50}, {42, -50}, {42, -32}}, color = {0, 127, 255}));
    connect(volume.ports[2], stopValve1.port_a) annotation(
      Line(points = {{42, -6}, {80, -6}}, color = {0, 127, 255}));
    connect(const.y, wellPump1.N_in) annotation(
      Line(points = {{-115, 56}, {-115, 55}, {-106, 55}, {-106, 20}}, color = {0, 0, 127}));
    connect(pipe1.port_b, tank.ports[1]) annotation(
      Line(points = {{-78, -30}, {-78, -50}, {-46, -50}, {-46, -30}}, color = {0, 127, 255}));
    connect(pipe2.port_b, volume1.ports[1]) annotation(
      Line(points = {{-136, -24}, {-136, 10}}, color = {0, 127, 255}));
    connect(boundary.ports[1], pipe2.port_a) annotation(
      Line(points = {{-150, -80}, {-136, -80}, {-136, -44}}, color = {0, 127, 255}));
    connect(wellPump1.port_b, pipe1.port_a) annotation(
      Line(points = {{-96, 10}, {-78, 10}, {-78, -10}}, color = {0, 127, 255}));
    connect(volume1.ports[2], wellPump1.port_a) annotation(
      Line(points = {{-136, 10}, {-116, 10}}, color = {0, 127, 255}));
    connect(ramp1.y, stopValve1.opening) annotation(
      Line(points = {{73, 32}, {90, 32}, {90, 2}}, color = {0, 0, 127}));
    connect(stopValve1.port_b, boundary1.ports[1]) annotation(
      Line(points = {{100, -6}, {120, -6}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 140, Tolerance = 1e-06, Interval = 0.28),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "ida", nls = "hybrid"),
      Diagram(coordinateSystem(extent = {{-200, -100}, {160, 100}})),
      Icon(coordinateSystem(extent = {{-200, -100}, {160, 100}})),
      __OpenModelica_commandLineOptions = "");
  end PressureTankTest2;

  model PressureSwitch
    replaceable package Medium = Modelica.Media.Water.StandardWater;
    ShallowWellPumpingSystem.StopValve stopValve1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {90, -6}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Sources.Boundary_pT boundary1(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {130, -6}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
    inner Modelica.Fluid.System system annotation(
      Placement(visible = true, transformation(origin = {-160, 46}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShallowWellPumpingSystem.PressureTank tank(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-46, -10}, extent = {{-20, -20}, {20, 20}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe(redeclare package Medium = Medium, diameter = 0.015, height_ab = 2, length = 2) annotation(
      Placement(visible = true, transformation(origin = {42, -22}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {42, 4}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    ShallowWellPumpingSystem.WellPump wellPump1(redeclare package Medium = Medium) annotation(
      Placement(visible = true, transformation(origin = {-106, 10}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe1(redeclare package Medium = Medium, diameter = 0.025, height_ab = -0.22, length = 0.22) annotation(
      Placement(visible = true, transformation(origin = {-78, -20}, extent = {{-10, -10}, {10, 10}}, rotation = -90)));
    Modelica.Fluid.Sources.Boundary_pT boundary(redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325) annotation(
      Placement(visible = true, transformation(origin = {-160, -80}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Fluid.Pipes.StaticPipe pipe2(redeclare package Medium = Medium, diameter = 0.06, height_ab = 6.3, length = 6.3) annotation(
      Placement(visible = true, transformation(origin = {-136, -34}, extent = {{-10, -10}, {10, 10}}, rotation = 90)));
    Modelica.Fluid.Vessels.ClosedVolume volume1(redeclare package Medium = Medium, V = 0.001, energyDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, massDynamics = Modelica.Fluid.Types.Dynamics.FixedInitial, nPorts = 2, p_start = system.p_start - 1000 * 6.3 * Modelica.Constants.g_n, use_portsData = false) annotation(
      Placement(visible = true, transformation(origin = {-136, 20}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Logical.Hysteresis hysteresis1(uHigh = 2.5e5, uLow = 8e4) annotation(
      Placement(visible = true, transformation(origin = {-36, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Logical.Not not1 annotation(
      Placement(visible = true, transformation(origin = {-4, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Logical.TriggeredTrapezoid triggeredTrapezoid1(amplitude = 3000, falling = 1, offset = 0, rising = 1) annotation(
      Placement(visible = true, transformation(origin = {28, 52}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
    Modelica.Blocks.Sources.Trapezoid trapezoid1(amplitude = 1, falling = 5, offset = 0, period = 120, rising = 5, startTime = 60, width = 30) annotation(
      Placement(visible = true, transformation(origin = {68, 34}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  equation
    connect(trapezoid1.y, stopValve1.opening) annotation(
      Line(points = {{80, 34}, {90, 34}, {90, 2}, {90, 2}}, color = {0, 0, 127}));
    connect(triggeredTrapezoid1.y, wellPump1.N_in) annotation(
      Line(points = {{40, 52}, {50, 52}, {50, 76}, {-106, 76}, {-106, 20}, {-106, 20}}, color = {0, 0, 127}));
    connect(not1.y, triggeredTrapezoid1.u) annotation(
      Line(points = {{8, 52}, {16, 52}}, color = {255, 0, 255}));
    connect(not1.u, hysteresis1.y) annotation(
      Line(points = {{-16, 52}, {-25, 52}}, color = {255, 0, 255}));
    connect(hysteresis1.u, tank.tankPressure) annotation(
      Line(points = {{-48, 52}, {-62, 52}, {-62, 12}}, color = {0, 0, 127}));
    connect(pipe.port_b, volume.ports[1]) annotation(
      Line(points = {{42, -12}, {42, -6}}, color = {0, 127, 255}));
    connect(tank.ports[2], pipe.port_a) annotation(
      Line(points = {{-46, -30}, {-40, -30}, {-40, -50}, {42, -50}, {42, -32}}, color = {0, 127, 255}));
    connect(volume.ports[2], stopValve1.port_a) annotation(
      Line(points = {{42, -6}, {80, -6}}, color = {0, 127, 255}));
    connect(pipe1.port_b, tank.ports[1]) annotation(
      Line(points = {{-78, -30}, {-78, -50}, {-46, -50}, {-46, -30}}, color = {0, 127, 255}));
    connect(pipe2.port_b, volume1.ports[1]) annotation(
      Line(points = {{-136, -24}, {-136, 10}}, color = {0, 127, 255}));
    connect(boundary.ports[1], pipe2.port_a) annotation(
      Line(points = {{-150, -80}, {-136, -80}, {-136, -44}}, color = {0, 127, 255}));
    connect(wellPump1.port_b, pipe1.port_a) annotation(
      Line(points = {{-96, 10}, {-78, 10}, {-78, -10}}, color = {0, 127, 255}));
    connect(volume1.ports[2], wellPump1.port_a) annotation(
      Line(points = {{-136, 10}, {-116, 10}}, color = {0, 127, 255}));
    connect(stopValve1.port_b, boundary1.ports[1]) annotation(
      Line(points = {{100, -6}, {120, -6}}, color = {0, 127, 255}));
    annotation(
      experiment(StartTime = 0, StopTime = 1000, Tolerance = 1e-06, Interval = 1),
      __OpenModelica_simulationFlags(lv = "LOG_STATS", outputFormat = "mat", s = "ida", nls = "hybrid"),
      Diagram(coordinateSystem(extent = {{-200, -100}, {160, 100}})),
      Icon(coordinateSystem(extent = {{-200, -100}, {160, 100}})),
      __OpenModelica_commandLineOptions = "");
  end PressureSwitch;
  annotation(
    uses(Modelica(version = "3.2.3")));
end ShallowWellPumpingSystem;
