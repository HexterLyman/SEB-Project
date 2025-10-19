import matplotlib.pyplot as plt
import re


with open("output.txt", "r") as file:
    cpp_output = file.read()


coordinates = re.findall(r"\((-?\d+\.?\d*),\s*(-?\d+\.?\d*)\)", cpp_output)
new_coordinate = re.search(r"New coordinates are: \((-?\d+\.?\d*),\s*(-?\d+\.?\d*)\)", cpp_output)
radius_match = re.search(r"radius is (\d+\.?\d*)", cpp_output)


coordinates = [(float(x), float(y)) for x, y in coordinates]
new_coordinate = (float(new_coordinate.group(1)), float(new_coordinate.group(2)))
radius = float(radius_match.group(1))


fig, ax = plt.subplots()
ax.scatter(*zip(*coordinates), label="Original Points", color="blue")
ax.scatter(*new_coordinate, color="red", label="Center Point")


circle = plt.Circle(new_coordinate, radius, color='b', fill=False, linestyle="--", label=f"Radius {radius}")
ax.add_patch(circle)


ax.set_aspect('equal', adjustable='box')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.xlabel("X-coordinate")
plt.ylabel("Y-coordinate")
plt.title("SEB_Project")
plt.grid(True)


plt.tight_layout()
plt.show()
