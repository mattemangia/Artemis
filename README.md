# Artemis Physics Engine

Un motore fisico 2D leggero scritto in C# con demo interattiva.

## Struttura del Progetto

### ArtemisEngine
Motore fisico 2D che include:
- Sistema di corpi rigidi (RigidBody)
- Forme fisiche (cerchi e rettangoli)
- Rilevamento e risoluzione delle collisioni
- Simulazione della gravità
- Fisica realistica con rimbalzo e attrito

### PhysicsCatapultDemo
Gioco demo che dimostra le capacità del motore Artemis:
- **Gameplay**: Lancia proiettili da una catapulta per distruggere strutture
- **Fisica Realistica**: Tutti gli oggetti seguono le leggi della fisica
- **Materiali Diversi**: Blocchi in legno, pietra, vetro e metallo con proprietà fisiche uniche
- **Sistema di Danno**: I blocchi si danneggiano in base all'impatto
- **Tre Livelli**: Strutture diverse da abbattere

### AdvancedPhysicsDemo
Demo interattivo che mostra le funzionalità avanzate del motore Artemis:

**5 Scene Interattive:**

1. **Ragdoll Physics** - Sistema di ragdoll con articolazioni
   - RevoluteJoint e DistanceJoint per braccia e gambe
   - Fisica realistica del corpo umano
   - Spawn multipli di ragdoll

2. **Vehicle with Suspension** - Veicolo con sospensioni
   - SpringJoint per ammortizzatori realistici
   - Controllo accelerazione e salto
   - Terreno irregolare con ostacoli

3. **Chains & Bridges** - Catene e ponti sospesi
   - DistanceJoint per catene metalliche
   - RopeJoint per ponti flessibili
   - Fisica realistica della tensione

4. **Raycast Shooting** - Sistema di raycasting
   - Visualizzazione raggi in tempo reale
   - Layer mask filtering
   - Applicazione impulsi ai punti di impatto

5. **Trigger Zones & Events** - Zone trigger ed eventi
   - Collision events (OnEnter, OnStay, OnExit)
   - Trigger zones senza fisica
   - Monitoraggio eventi in tempo reale

## Controlli

- **W/S**: Regola l'angolo di lancio (su/giù)
- **A/D**: Regola la potenza del lancio (meno/più)
- **SPAZIO**: Lancia un proiettile
- **R**: Ricarica il livello corrente
- **1/2/3**: Carica il livello 1/2/3
- **Q**: Esci dal gioco

## Come Compilare ed Eseguire

```bash
# Compila l'intera solution
dotnet build Artemis.sln

# Esegui Physics Catapult Demo
dotnet run --project PhysicsCatapultDemo/PhysicsCatapultDemo.csproj

# Esegui Advanced Physics Demo
dotnet run --project AdvancedPhysicsDemo/AdvancedPhysicsDemo.csproj
```

## Caratteristiche del Motore Artemis

### Sistema Fisico Core
- Integrazione di Eulero per il movimento
- Rilevamento collisioni tra cerchi, rettangoli e forme miste
- Risoluzione delle collisioni con impulsi
- Attrito e rimbalzo configurabili
- Correzione della penetrazione
- Momento di inerzia per rotazioni realistiche
- Applicazione di forze, impulsi e torque

### Raycasting
- Sistema completo di raycast per cerchi e rettangoli
- Supporto per layer mask
- RaycastAll per rilevare tutti gli oggetti su un raggio
- Calcolo di punto di impatto, normale e distanza

### Sistema di Collisione Avanzato
- **Collision Layers**: Filtraggio bit-mask per controllare cosa collide con cosa
- **Collision Events**: Eventi OnCollisionEnter, OnCollisionStay, OnCollisionExit
- **Trigger Zones**: Sensori che rilevano collisioni senza fisica
- **Trigger Events**: OnTriggerEnter, OnTriggerStay, OnTriggerExit
- **Collision Listeners**: Interface per callback personalizzati

### Ottimizzazione Performance
- **Spatial Hash Grid**: Partizione spaziale per collisioni O(n) invece di O(n²)
- **QuadTree**: Alternativa al grid per distribuzioni non uniformi
- **Sleeping System**: Corpi inattivi vengono messi in "sleep" per risparmiare CPU
  - Wake up automatico su impulsi o collisioni
  - Soglie configurabili di velocità e tempo

### Joints e Constraints
- **DistanceJoint**: Mantiene distanza fissa tra due corpi
- **SpringJoint**: Molla con costante elastica e smorzamento
- **RevoluteJoint**: Rotazione attorno a un punto fisso (cerniera)
- **RopeJoint**: Vincolo di distanza massima
- Tutti i joint supportano stiffness configurabile

### Continuous Collision Detection (CCD)
- **Swept Collision**: Previene il tunneling per oggetti veloci
- **Time of Impact (TOI)**: Calcolo preciso del momento d'impatto
- **Conservative Advancement**: Algoritmo iterativo per forme complesse
- Abilitabile per singoli corpi con `UseCCD = true`

### Forme Avanzate
- **PolygonShape**: Poligoni convessi arbitrari (triangoli, pentagoni, esagoni, ecc.)
- **CompoundShape**: Forme composite da più shape semplici
- **EdgeShape**: Segmenti per terreni statici
- Creazione procedurale di poligoni regolari

### Area Effectors
Effetti ambientali che influenzano i corpi in un'area:
- **DirectionalForce**: Vento, correnti d'acqua, gravità direzionale
- **RadialForce**: Pozzi gravitazionali, esplosioni, campi magnetici
- **Buoyancy**: Simulazione acqua con galleggiamento e resistenza
- **Vortex**: Tornado, mulinelli, vortici
- **Point**: Forze puntuali (attrazioni/repulsioni)
- Filtri per layer e tag per controllo preciso

### One-Way Platforms
- Piattaforme attraversabili da una direzione
- **JumpThrough**: Piattaforme su cui si può saltare dal basso
- **FallThrough**: Piattaforme che si possono attraversare dall'alto
- Configurabili con soglia di velocità personalizzabile

### Solver Avanzato
- **Sequential Impulse Solver**: Risoluzione iterativa per maggiore stabilità
- **Contact Persistence**: Mantiene i contatti tra frame
- **Warm Starting**: Usa impulsi precedenti per convergenza più veloce
- **Coulomb Friction**: Modello di attrito realistico con cono di attrito
- **Baumgarte Stabilization**: Correzione della penetrazione graduale
- Configurabile: numero di iterazioni velocità/posizione, tolleranza

### Fisica Scientifica
Sistema completo per applicazioni scientifiche e simulazioni ad alta precisione:

#### Metodi di Integrazione Avanzati
- **Euler**: Integrazione base (veloce ma meno precisa)
- **Semi-Implicit Euler**: Stabile ed efficiente per giochi
- **Verlet**: Migliore conservazione dell'energia
- **Runge-Kutta 4 (RK4)**: Massima precisione per simulazioni scientifiche

#### Energy Tracking
- Monitoraggio energia cinetica, potenziale e rotazionale
- Calcolo drift energetico (perdita/guadagno energia nel tempo)
- Percentuale di deriva per validare la simulazione

#### Deterministic Physics
- Simulazioni riproducibili con seed fisso
- Stesse condizioni iniziali = stessi risultati
- Essenziale per testing e replay

#### Data Export
- Esportazione dati in CSV per analisi
- Registrazione posizione, velocità, rotazione, energia per frame
- Integrabile con Python, MATLAB, Excel per post-processing

### Materiali
Ogni materiale ha proprietà uniche:

- **Legno**: Leggero, fragile, alta attrito
- **Pietra**: Pesante, resistente, molto attrito
- **Vetro**: Pesante ma fragile, basso attrito
- **Metallo**: Molto pesante, molto resistente

### Sistema di Rendering
Rendering ASCII/testuale per console che visualizza:
- Corpi fisici con caratteri diversi per materiale
- Indicatori di salute con cambio colore
- UI con statistiche di gioco
- Mirino della catapulta

## Architettura

Il progetto è diviso in due componenti principali:

1. **ArtemisEngine**: Libreria core del motore fisico
   - `Vector2.cs`: Matematica vettoriale 2D
   - `RigidBody.cs`: Corpi rigidi e forme con proprietà fisiche avanzate
   - `PhysicsWorld.cs`: Simulazione del mondo fisico con spatial partitioning
   - `CollisionDetection.cs`: SAT e algoritmi di collision detection avanzati
   - `Raycasting.cs`: Sistema di raycast per linee di vista e proiettili
   - `CollisionLayers.cs`: Sistema di layer per filtraggio collisioni
   - `CollisionEvents.cs`: Eventi e callback per collisioni
   - `SleepingSystem.cs`: Ottimizzazione per corpi inattivi
   - `Joints.cs`: Constraints e collegamenti tra corpi
   - `SpatialPartitioning.cs`: Spatial hash grid e quadtree
   - `ContinuousCollision.cs`: CCD per prevenire tunneling
   - `AdvancedShapes.cs`: Poligoni, forme composite, edge shapes
   - `AreaEffectors.cs`: Effetti ambientali (vento, gravità, buoyancy)
   - `OneWayPlatform.cs`: Piattaforme attraversabili
   - `ConstraintSolver.cs`: Sequential Impulse solver con warm starting
   - `ScientificPhysics.cs`: Integratori avanzati, energy tracking, data export

2. **PhysicsCatapultDemo**: Gioco dimostrativo
   - `Material.cs`: Definizione dei materiali
   - `GameObject.cs`: Oggetti di gioco con salute
   - `Catapult.cs`: Sistema di lancio
   - `LevelBuilder.cs`: Costruttore di livelli
   - `Renderer.cs`: Sistema di rendering ASCII
   - `Program.cs`: Game loop principale

3. **AdvancedPhysicsDemo**: Demo funzionalità avanzate
   - `RagdollBuilder.cs`: Creazione ragdoll con joints
   - `VehicleBuilder.cs`: Veicolo con sospensioni
   - `ChainBuilder.cs`: Catene e ponti sospesi
   - `AdvancedRenderer.cs`: Rendering con visualizzazione joints e raycast
   - `Program.cs`: 5 scene interattive

## Esempi di Utilizzo

### Creare un Corpo Rigido
```csharp
// Creare una sfera
var circleShape = new CircleShape(1.0f);
var ball = new RigidBody(new Vector2(10, 10), mass: 5.0f, circleShape);
ball.Restitution = 0.8f; // Rimbalzo
ball.Friction = 0.3f;

// Creare un rettangolo
var boxShape = new BoxShape(width: 2.0f, height: 1.0f);
var box = new RigidBody(new Vector2(5, 5), mass: 10.0f, boxShape);

// Aggiungere al mondo
world.AddBody(ball);
world.AddBody(box);
```

### Collision Layers
```csharp
// Impostare layer e mask
projectile.CollisionLayer = CollisionLayers.Projectile;
projectile.CollisionMask = CollisionLayers.Enemy | CollisionLayers.Ground;

enemy.CollisionLayer = CollisionLayers.Enemy;
enemy.CollisionMask = CollisionLayers.Projectile | CollisionLayers.Player;
```

### Raycasting
```csharp
var ray = new Ray(origin: new Vector2(0, 0),
                  direction: new Vector2(1, 0),
                  maxDistance: 100f);

var hit = world.Raycast(ray, layerMask: CollisionLayers.Enemy);
if (hit.Hit)
{
    Console.WriteLine($"Hit at {hit.Point}, distance: {hit.Distance}");
}
```

### Joints
```csharp
// Collegare due corpi con una molla
var spring = new SpringJoint(bodyA, bodyB,
                             restLength: 5.0f,
                             springConstant: 100f,
                             damping: 0.5f);
world.AddJoint(spring);

// Vincolo di distanza
var distance = new DistanceJoint(bodyA, bodyB, distance: 3.0f);
world.AddJoint(distance);
```

### Trigger Zones
```csharp
// Creare una zona trigger
var trigger = new RigidBody(new Vector2(20, 20), 0, new CircleShape(5f));
trigger.IsTrigger = true; // Non ha fisica, solo detection
trigger.IsStatic = true;

// Ascoltare eventi
world.OnTriggerEnter += (sender, e) =>
{
    Console.WriteLine($"Oggetto entrato nella trigger zone!");
};
```

### Collision Events
```csharp
world.OnCollisionEnter += (sender, e) =>
{
    Console.WriteLine($"Collisione tra {e.BodyA} e {e.BodyB}");
    // Applicare danno, suoni, effetti, ecc.
};
```

### Continuous Collision Detection (CCD)
```csharp
// Abilitare CCD per proiettili veloci (previene tunneling)
var bullet = new RigidBody(new Vector2(0, 0), 0.1f, new CircleShape(0.2f));
bullet.UseCCD = true; // Usa swept collision detection
bullet.Velocity = new Vector2(100, 0); // Velocità molto alta

world.AddBody(bullet);
```

### Forme Avanzate
```csharp
// Creare un poligono regolare (esagono)
var hexagon = PolygonShape.CreateRegular(sides: 6, radius: 2.0f);
var hexBody = new RigidBody(new Vector2(10, 10), 5f, hexagon);

// Creare una forma composita (es. corpo umano)
var compound = new CompoundShape();
compound.AddShape(new CircleShape(1f), new Vector2(0, 2)); // Testa
compound.AddShape(new BoxShape(1.5f, 2f), new Vector2(0, 0)); // Torso
var character = new RigidBody(new Vector2(5, 5), 10f, compound);

world.AddBody(hexBody);
world.AddBody(character);
```

### Area Effectors
```csharp
// Creare una zona con vento
var wind = new DirectionalForceEffector(
    center: new Vector2(20, 20),
    radius: 10f,
    force: new Vector2(50, 0), // Vento verso destra
    drag: 0.1f
);
world.AddAreaEffector(wind);

// Creare una zona di galleggiamento (acqua)
var water = new BuoyancyEffector(
    center: new Vector2(30, 5),
    radius: 15f,
    density: 1.0f, // Densità dell'acqua
    linearDrag: 2.0f,
    angularDrag: 0.5f
);
water.FlowVelocity = new Vector2(2, 0); // Corrente
world.AddAreaEffector(water);

// Pozzo gravitazionale (buco nero)
var gravityWell = new RadialForceEffector(
    center: new Vector2(50, 50),
    radius: 20f,
    strength: -100f, // Negativo = attrazione
    falloff: RadialFalloff.InverseSquare
);
world.AddAreaEffector(gravityWell);
```

### One-Way Platforms
```csharp
// Creare piattaforma jump-through (salta dal basso)
var platform = new RigidBody(new Vector2(20, 10), 0, new BoxShape(10f, 1f));
platform.IsStatic = true;
platform.SetOneWayPlatform(OneWayPlatformBehavior.CreateJumpThrough());

world.AddBody(platform);
```

### Solver Avanzato e Fisica Scientifica
```csharp
// Configurare il solver avanzato
world.UseAdvancedSolver = true; // Abilita Sequential Impulse solver

// Usare integrazione RK4 per massima precisione
var body = new RigidBody(new Vector2(0, 50), 10f, new CircleShape(1f));
world.AddBody(body);

// Definire funzione accelerazione personalizzata
Func<Vector2, Vector2, Vector2> accelFunc = (pos, vel) =>
{
    // Gravità + attrito aria
    Vector2 gravity = new Vector2(0, -9.81f);
    Vector2 drag = -vel * 0.1f;
    return gravity + drag;
};

// Applicare RK4 nel game loop
AdvancedIntegration.RK4Integration(body, accelFunc, deltaTime);

// Monitorare energia per validare conservazione
var energyTracker = new EnergyTracker();
energyTracker.GravityMagnitude = 9.81f;
energyTracker.UpdateEnergy(world.Bodies);

float initialEnergy = energyTracker.TotalEnergy;
// ... simulazione ...
energyTracker.UpdateEnergy(world.Bodies);
float energyDrift = energyTracker.GetEnergyDriftPercentage(initialEnergy);
Console.WriteLine($"Energy drift: {energyDrift:F2}%");

// Esportare dati per analisi
var exporter = new SimulationDataExporter();
// Nel game loop:
exporter.RecordFrame(world.Bodies, deltaTime);
// Alla fine:
string csv = exporter.ExportToCSV();
File.WriteAllText("simulation_data.csv", csv);
```

## Requisiti

- .NET 8.0 SDK o superiore
- Console con supporto per colori ANSI

## Licenza

Progetto dimostrativo per il motore fisico Artemis.
