import java.io.PrintWriter;
import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

public class Simulation {

    private static final int    RUN = 10;

    private static final double DESIRED_VEL = 2;
    private static final boolean THICC_MODE = false;
    private static final int    BASE = 1;                       // DT base
    private static final int    EXP = 5;                        // DT exp
    private static final double DT = BASE * Math.pow(10, -EXP); // Step delta time
    private static final int    N = 200;                        // Number of particles
    private static final double WIDTH = 20;
    private static final double HEIGHT = 20;
    private static final double SLIT_SIZE = 1.2;
    private static final double k = 1.2e5;
    private static final double kt = 2.4e5;
    private static final double gamma = 140;
    private static final double A = 2000;
    private static final double B = 0.08;
    private static final double tau = 0.5;
    private static final double SLIT_Xo = (WIDTH - SLIT_SIZE)/2;
    private static final double SLIT_Xf = ((WIDTH - SLIT_SIZE)/2) + SLIT_SIZE;
    private static final double MIN_PARTICLE_R = 0.25;          // Min particle radius
    private static final double MAX_PARTICLE_R = 0.29;         // Max particle radius
    private static final double STEP_PRINT_DT = 0.1;
    private static final double ANIMATION_DT = 1.0 / 60;          // DT to save a simulation state
    private static final double MEASURE_DT = 1.0 / 10;                // DT to save a simulation state
    private static final double MAX_SIM_TIME = 10;             // Max simulation time in seconds

    private static double              simTime = 0; //Simulation time in seconds
    private static List<Particle> particles = new ArrayList<>(N);
    private static ArrayList<Wall>     walls = new ArrayList<>(4);

    private static List<List<Particle>> savedStates = new ArrayList<>();
    private static List<Double> kineticEnergy = new LinkedList<>();
    private static List<Double> times = new LinkedList<>();

    public static void main(String[] args) throws Exception{
        System.out.println(String.format("Run ID: %d - N: %d", RUN, N));
        PrintWriter writer = new PrintWriter("data/" + DESIRED_VEL + "_" + gamma + "_" + BASE + "e-" + EXP + "_" + THICC_MODE + "_simulation_" + RUN + ".xyz");

        initWalls(WIDTH, HEIGHT, SLIT_SIZE);
        initParticles(N, WIDTH, HEIGHT, MIN_PARTICLE_R, MAX_PARTICLE_R);

        saveMeasures();
        writeState(writer);

        int lastFrame = 1, lastMeasure = 1, lastStepPrint = 0;
        System.out.println("Starting simulation");

        List<Double> exitTimes = new LinkedList<>();

        while(particles.size() > 0) {
            // Clear forces and add interaction forces with walls to particles and add driving force
            particles.forEach(p -> {
                p.clearForces();
                applyDrivingForce(p);
                for (Wall w : walls) {
                    double overlap = w.getOverlap(p);
                    if (overlap > 0) {
                        applyGranularForce(w, p, overlap);
                    }
                }
            });

            // Add interaction forces between particles
            IntStream.range(0, particles.size()).forEach(i -> {
                Particle pi = particles.get(i);

                for (int j = i + 1; j < particles.size(); j++) {
                    Particle pj = particles.get(j);
                    double overlap = pi.getOverlap(pj);
                    if (overlap > 0) {
                        applyGranularForce(pi, pj, overlap);
                    }
                    applySocialForce(pi, pj);
                }
            });

            // Move particles a DT time and filter the ones that are out
            particles = particles.stream().peek(p -> {
                p.move(DT);
            }).filter(Simulation::isInScene).collect(Collectors.toList());

            // Add DT to simulation time
            simTime += DT;

            particles.stream().filter(Simulation::isOut).forEach((p) -> {
                    exitTimes.add(simTime);
                    p.isOut = true;
                });


            if (simTime / STEP_PRINT_DT > lastStepPrint) {
                System.out.println(String.format("simTime: %.2f", simTime));
                lastStepPrint++;
            }

            if (simTime / ANIMATION_DT > lastFrame) {
                writeState(writer);
                lastFrame++;
            }

            if (simTime / MEASURE_DT > lastMeasure) {
                saveMeasures();
                lastMeasure++;
            }
        }
        saveMeasures();
        System.out.println("Finished simulation");
        writer.close();

        System.out.println("Printing measures");
        System.out.println(String.format("Evacuated particles: %d", exitTimes.size()));

        printList(kineticEnergy, "data/" + DESIRED_VEL + "_" + gamma + "_" + BASE + "e-" + EXP + "_" + THICC_MODE + "_kineticEnergy_" + RUN + ".csv");
        printList(times, "data/" + DESIRED_VEL + "_" + gamma + "_" + BASE + "e-" + EXP + "_" + THICC_MODE + "_times_" + RUN + ".csv");
        printList(exitTimes, "data/" + DESIRED_VEL + "_" + gamma + "_" + BASE + "e-" + EXP + "_" + THICC_MODE + "_exitTimes_" + RUN + ".csv");

    }

    private static boolean isOut(Particle p) {
        return !p.isOut && p.y < 0;
    }

    private static boolean isInScene(Particle p) {
        return p.y > - HEIGHT / 10;
    }

    private static void applyDrivingForce(Particle p) {
        double targetX;
        double targetY = 0;
        double minTargetX =  SLIT_Xo + p.r;
        double maxTargetX =  SLIT_Xf - p.r;
        if(p.x >= minTargetX && p.x <= maxTargetX){
            targetX = p.x;
        }else if(p.x > maxTargetX){
            targetX = maxTargetX;
        }else{
            targetX = minTargetX;
        }
        if(p.y <= 0){
            targetY = - HEIGHT / 10;
        }

        double eTargetX = p.eTargetX(targetX, targetY);
        double eTargetY = p.eTargetY(targetX, targetY);

        double fx = p.m * (DESIRED_VEL * eTargetX - p.vx)/tau;
        double fy = p.m * (DESIRED_VEL * eTargetY - p.vy)/tau;

        p.fx += fx;
        p.fy += fy;
    }

    private static void applySocialForce(Particle p1, Particle p2) {

        double enx = p1.enx(p2);
        double eny = p1.eny(p2);

        double distanceBetweenBorders = p1.centerDistance(p2) - p1.r - p2.r;

        double fn = -A*Math.exp(-distanceBetweenBorders/B);

        double fx = fn * enx;
        double fy = fn * eny;

        p1.fx += fx;
        p1.fy += fy;
        p2.fx -= fx;
        p2.fy -= fy;
    }

    private static void applyGranularForce(Particle p1, Particle p2, double overlap) {

        double enx = p1.enx(p2);
        double eny = p1.eny(p2);

        double tangentRelVel = p1.getTangentRelVel(p2);

        double ft = -kt*overlap*tangentRelVel;
        double fn = -k*overlap;

        double fx = fn * enx -eny*ft;
        double fy = fn * eny + enx*ft;

        double fn_mod = Math.abs(fn);
        p1.fn += fn_mod;
        p1.fx += fx;
        p1.fy += fy;
        p2.fn += fn_mod;
        p2.fx -= fx;
        p2.fy -= fy;
    }

    private static void applyGranularForce(Wall w, Particle p, double overlap) {

        double normalRelVel = w.getNormalRelVel(p);

        double fn = -k*overlap - gamma*normalRelVel;

        double fx = fn * w.enx(p);
        double fy = fn * w.eny(p);

        p.fn += Math.abs(fn);
        p.fx += fx;
        p.fy += fy;
    }

    private static void initWalls(double width, double height, double slitSize) {
        walls.add(new Wall(0, 0, 0, height, -1, 0));
        walls.add(new Wall(width, 0, width, height, 1, 0));
        walls.add(new Wall(0, height, width, height, 0, 1));
        if(slitSize == 0) {
            walls.add(new Wall(0, 0, width, 0, 0, -1));
        }else{
            double bottomWallWidth = (width - slitSize) / 2;
            walls.add(new Wall(0, 0, bottomWallWidth, 0, 0, -1));
            walls.add(new Wall(width - bottomWallWidth, 0, width, 0, 0, -1));
        }
    }

    private static void initParticles(int n, double width, double height, double minRadius, double maxRadius) {
        while (particles.size() < n) {
            double particleRadius = minRadius + Math.random() * (maxRadius - minRadius);
            double x = particleRadius + Math.random() * (width - 2 * particleRadius);
            double y = particleRadius + Math.random() * (height - 2 * particleRadius);

            Particle newParticle = new Particle(particles.size(), x, y, particleRadius, THICC_MODE);

            boolean valid = particles.stream().parallel().allMatch(p -> p.getOverlap(newParticle) == 0);

            if (valid) {
                particles.add(newParticle);
            }
        }
    }

    private static void printList(List<Double> list, String filename) {
        try {
            PrintWriter writer = new PrintWriter(filename);
            list.forEach(writer::println);
            writer.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    private static void saveMeasures() {
        times.add(simTime);
        kineticEnergy.add(particles.parallelStream().map(Particle::kineticEnergy).reduce(0.0, (d1, d2) -> d1 + d2));
    }

    private static void writeState(PrintWriter writer) {
        writer.println(particles.size() + 2);
        writer.println();
        writer.println("-2 0.0 0.0 0.00000001 0.0 0.0");
        writer.println(String.format(Locale.ENGLISH, "-1 %f %f 0.00000001 0.0 0.0", WIDTH, HEIGHT));
        particles.stream().parallel().forEach(writer::println);
    }
}
