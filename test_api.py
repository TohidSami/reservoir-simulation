import requests

url = "http://127.0.0.1:5000/api/simulation"

# داده‌های ورودی (JSON)
payload = {
    "name": "API Test Case 1",
    "nx": 15, "ny": 15, "nz": 3,
    "swc": 0.12, "sor": 0.25,
    "inj_rate_target": 2000.0,
    "prod_bhp_min": 2500.0,
    "dt": 10,
    "n_steps": 10
}

print("Sending request to server...")
response = requests.post(url, json=payload)

print(f"Status Code: {response.status_code}")
print("Response:", response.json())