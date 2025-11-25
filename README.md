# 🛢️ Reservoir Simulation REST API

A Dockerized RESTful API for numerical reservoir simulation, built with Python, Flask, and PostgreSQL.

## 🚀 Overview
This project transforms a traditional reservoir engineering script into a scalable web service. It allows users to:
1.  Submit simulation parameters (Grid, PVT, Wells) via API.
2.  Run a fully implicit numerical solver (Newton-Raphson) in the background.
3.  Store results in a structured PostgreSQL database.
4.  Retrieve production plots and pressure maps via API endpoints.

## 🛠️ Tech Stack
* **Language:** Python 3.9
* **Framework:** Flask (Backend API)
* **Database:** PostgreSQL (via SQLAlchemy ORM)
* **Simulation Engine:** Custom NumPy-based solver (3D, 2-Phase Flow)
* **Infrastructure:** Docker & Docker Compose
* **Visualization:** Matplotlib

## ⚙️ How to Run
Prerequisites: Docker Desktop installed.

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/YOUR_USERNAME/reservoir-api.git](https://github.com/YOUR_USERNAME/reservoir-api.git)
    cd reservoir-api
    ```

2.  **Start Services:**
    ```bash
    docker-compose up --build
    ```

3.  **Access API:**
    The server will start at `http://localhost:5000`.

## 🔌 API Endpoints
| Method | Endpoint | Description |
| :--- | :--- | :--- |
| `POST` | `/api/simulation` | Run a new simulation case |
| `GET` | `/api/simulation/<id>` | Get simulation results (Pressure/Saturation matrices) |
| `GET` | `/api/simulation/<id>/well_plot` | Generate Production vs. Injection plot |

## 📊 Architecture
The project follows a **Containerized Architecture**:
* **Service A (Web):** Flask Application handling logic and simulation engine.
* **Service B (DB):** PostgreSQL container for persistent storage.
* **Networking:** Docker Compose internal network for secure communication.