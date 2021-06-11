
model ValvesLinearAndOnOff
  replaceable package Medium = Modelica.Media.Water.StandardWater;
  Modelica.Blocks.Sources.Sine sine(amplitude = 0.5, freqHz = 0.5, offset = 0.5, phase = 0, startTime = 0)  annotation(
    Placement(visible = true, transformation(origin = {-32, 66}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT boundary(
  redeclare package Medium = Medium, T = 293.15, nPorts = 2, p = 121325)  annotation(
    Placement(visible = true, transformation(origin = {-70, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT boundary2(
  redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325)  annotation(
    Placement(visible = true, transformation(origin = {68, -30}, extent = {{10, -10}, {-10, 10}}, rotation = 0)));
  Modelica.Fluid.Sources.Boundary_pT boundary_pT(
  redeclare package Medium = Medium, T = 293.15, nPorts = 1, p = 101325)  annotation(
    Placement(visible = true, transformation(origin = {64, 24}, extent = {{10, 10}, {-10, -10}}, rotation = 0)));
  Modelica.Fluid.Valves.ValveLinear valveLinear(
  redeclare package Medium = Medium, dp_nominal = 10000, m_flow_nominal = 1)  annotation(
    Placement(visible = true, transformation(origin = {-2, 26}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Fluid.Valves.ValveDiscrete valveDiscrete(
  redeclare package Medium = Medium, dp_nominal = 10000, m_flow_nominal = 1, opening_min = 0.001)  annotation(
    Placement(visible = true, transformation(origin = {0, -30}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  Modelica.Blocks.Sources.BooleanPulse booleanPulse(period = 2, startTime = 0, width = 50)  annotation(
    Placement(visible = true, transformation(origin = {-30, -2}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
  inner Modelica.Fluid.System system annotation(
    Placement(visible = true, transformation(origin = {86, 86}, extent = {{-10, -10}, {10, 10}}, rotation = 0)));
equation
  connect(sine.y, valveLinear.opening) annotation(
    Line(points = {{-20, 66}, {-2, 66}, {-2, 34}, {-2, 34}}, color = {0, 0, 127}));
  connect(valveLinear.port_b, boundary_pT.ports[1]) annotation(
    Line(points = {{8, 26}, {54, 26}, {54, 24}, {54, 24}}, color = {0, 127, 255}));
  connect(valveLinear.port_a, boundary.ports[1]) annotation(
    Line(points = {{-12, 26}, {-60, 26}, {-60, 26}, {-60, 26}}, color = {0, 127, 255}));
  connect(boundary.ports[2], valveDiscrete.port_a) annotation(
    Line(points = {{-60, 26}, {-48, 26}, {-48, -30}, {-10, -30}, {-10, -30}}, color = {0, 127, 255}));
  connect(valveDiscrete.port_b, boundary2.ports[1]) annotation(
    Line(points = {{10, -30}, {58, -30}}, color = {0, 127, 255}));
  connect(booleanPulse.y, valveDiscrete.open) annotation(
    Line(points = {{-18, -2}, {0, -2}, {0, -22}, {0, -22}}, color = {255, 0, 255}));
end ValvesLinearAndOnOff;
